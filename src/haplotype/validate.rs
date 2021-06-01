// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
use std::path::Path;
use std::collections::HashSet;
use std::fs::File;
use std::cmp::max;
use std::cmp::min;
use std::str;
use std::str::FromStr;
use std::cmp::PartialOrd;
use std::ops::Range;
use std;

use serde::Serialize;
use csv;
use itertools::Itertools;

use bio::io::fasta;
use bio::alphabets;
use bio::alignment::pairwise::banded;
use bio::alignment::{self, Alignment, sparse};

use bio::stats::probs::LogProb;
use bio::stats::combinatorics::combinations;
use rust_htslib::bam::{self, Read, IndexedReader, Record};

use event::BedpeRow;
use locus::Locus;
use call;
use event;
use Args;

type Aligner = banded::Aligner<alignment::pairwise::MatchParams>;

fn get_read(rec: &Record) -> Vec<u8> {
    let raw_seq = rec.seq().as_bytes();

    if rec.is_reverse() {
        event::rc_seq(&raw_seq)
    } else {
        raw_seq
    }
}


fn aln_read(read_seq: &[u8],
            ref_seq: &[u8],
            rc_ref_seq: &[u8],
            aligner: &mut Aligner)
            -> (i32, usize, bool) {
    let fwd = aligner.custom(read_seq, ref_seq);
    let rev = aligner.custom(read_seq, rc_ref_seq);

    if fwd.score > rev.score {
        (fwd.score, fwd.ystart, true)
    } else {
        (rev.score, ref_seq.len() - rev.ystart, false)
    }
}

// Score for a read pair
fn score_read_pair(r1: &[u8], r2: &[u8], hap: &[u8], rc_hap: &[u8], aligner: &mut Aligner) -> i32 {

    let v1 = aln_read(r1, hap, rc_hap, aligner);
    let v2 = aln_read(r2, hap, rc_hap, aligner);
    let (score1, pos1, strand1) = v1;
    let (score2, pos2, strand2) = v2;

    let proper_pair = (strand1 != strand2) &&
        ((pos1 < pos2 && strand1) || (pos2 < pos1 && strand2));
    let proper_pair_score = if proper_pair { 0 } else { -10 };

    let dist = (pos1 as f64) - (pos2 as f64);
    let mut dist_score = -(dist - 300.0).powf(2.0) / (2.0 * 100.0f64.powf(2.0));
    if dist_score < -10.0 {
        dist_score = -10.0;
    }

    //println!("r1: {:?}, r2: {:?}", v1, v2);
    //println!("s1: {}, s2: {}, pp: {}, d:{}", score1, score2, proper_pair, dist_score);
    //let dist_score = 0;

    score1 + score2 + proper_pair_score + (dist_score as i32)
}

// Score for a single read
fn score_single_read(r1: &[u8], hap: &[u8], rc_hap: &[u8], aligner: &mut Aligner) -> i32 {
    let (score1, _, _) = aln_read(r1, hap, rc_hap, aligner);
    score1
}

fn score_hap_read_pair(rp: &HapReadPair, hap: &[u8], rc_hap: &[u8], aligner: &mut Aligner) -> i32 {
    match rp.seq2 {
        Some(ref s2) => score_read_pair(&rp.seq1, s2, hap, rc_hap, aligner),
        _ => score_single_read(&rp.seq1, hap, rc_hap, aligner),
    }
}


/// return Vector of sequences (seq, qual, bc, read)
fn get_reads(bam: &mut IndexedReader, locus: &Locus) -> Vec<Vec<u8>> {
    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();

    let mut seqs = Vec::new();

    // Add phased reads to haplotype buckets
    bam.fetch(tid, locus.start, locus.end)
        .expect("Error fetching BAM file.");
    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();

        // Skip poor mappings & secondary alignments
        if rec.mapq() < 10 || rec.is_secondary() {
            continue;
        }

        let raw_seq = call::get_seq_bounded(&rec, locus);

        let fwd_seq = if rec.is_reverse() {
            event::rc_seq(&raw_seq)
        } else {
            raw_seq
        };

        seqs.push(fwd_seq);
    }

    seqs
}

struct HapReadPair {
    hap: Option<u8>,
    hap_qual: Option<u8>,
    seq1: Vec<u8>,
    seq2: Option<Vec<u8>>,
}

/// return Vector of sequences (seq, qual, bc, read)
fn get_read_pairs(bam: &mut IndexedReader, locus: &Locus) -> Vec<HapReadPair> {
    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();

    // Add phased reads to haplotype buckets
    bam.fetch(tid, locus.start, locus.end)
        .expect("Error fetching BAM file.");
    let mut recs = vec![];

    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();
        // Skip poor mappings & secondary alignments
        if rec.mapq() <= 10 || rec.is_secondary() {
            continue;
        }
        recs.push(rec);
    }

    recs.sort_by(|x, y| x.qname().cmp(y.qname()));

    let mut seqs = Vec::new();
    for (_, pair) in &recs.iter().group_by(|x| x.qname()) {
        let mut pair_iter = pair.into_iter();
        let r1 = pair_iter.next();
        let r2 = pair_iter.next();

        let hap = call::get_haplotype(r1.unwrap());
        let hap_qual = call::get_haplotype_qual(r1.unwrap());

        let seq1 = r1.map(|x| {
            let mut s = x.seq().as_bytes();
            if x.is_reverse() {
                s = alphabets::dna::revcomp(&s);
            }

            s
        })
            .unwrap();

        let seq2 = r2.map(|x| {
            let mut s = x.seq().as_bytes();
            if x.is_reverse() {
                s = alphabets::dna::revcomp(&s);
            }

            s
        });
        let hrp = HapReadPair {
            hap,
            hap_qual,
            seq1,
            seq2,
        };
        seqs.push(hrp);
    }

    seqs
}

/// Validate the event using the coverage information. The algorithm works as follows
/// 1. Compute mean coverage around the locus
/// 2. Compute phased coverage within the event
/// 3. Score each hypothesis based on the phased coverage in the event region
pub fn validate_using_coverage
(bam: &mut IndexedReader,
 locus: &Locus,
 event: &(Option<((usize, usize), usize)>, Option<((usize, usize), usize)>))
 -> (Option<((usize, usize), usize)>, Option<((usize, usize), usize)>) {

    let min_mean_coverage = 5; // If the mean coverage is less than this, bail out

    let empty_range = Range {
        start: std::usize::MAX,
        end: 0,
    };

    // Each hypothesis is represented as a deletion range
    let (del_ranges, hypotheses) = match *event {
        (None, None) => panic!("don't validate empty events"),
        (None, Some(e)) => {
            let ((ev_len, ev_start), _) = e;
            let del_range = Range {
                start: ev_start,
                end: ev_start + ev_len,
            };
            (vec![empty_range, del_range], vec![(0, 0), (0, 1), (1, 1)])
        }
        (Some(e), None) => {
            let ((ev_len, ev_start), _) = e;
            let del_range = Range {
                start: ev_start,
                end: ev_start + ev_len,
            };
            (vec![empty_range, del_range], vec![(0, 0), (1, 0), (1, 1)])
        }
        (Some(e1), Some(e2)) if e1 == e2 => {
            let ((ev_len, ev_start), _) = e1;
            let del_range = Range {
                start: ev_start,
                end: ev_start + ev_len,
            };
            (vec![empty_range, del_range], vec![(0, 0), (0, 1), (1, 0), (1, 1)])
        }
        (Some(e1), Some(e2)) => {
            let ((ev_len, ev_start), _) = e1;
            let del_range1 = Range {
                start: ev_start,
                end: ev_start + ev_len,
            };
            let ((ev_len, ev_start), _) = e2;
            let del_range2 = Range {
                start: ev_start,
                end: ev_start + ev_len,
            };
            (vec![empty_range, del_range1, del_range2], vec![(0, 0), (1, 0), (0, 2), (1, 2)])
        }
    };

    // Compute the mean coverage

    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    // First check the mean read length within the locus to
    // get rid of boundary bias when doing pileup
    bam.fetch(tid, locus.start, locus.end)
        .expect("Error fetching BAM file.");
    let mut max_read_len = 0;
    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();
        max_read_len = max(max_read_len, rec.seq().len());
    }

    // how many bases to scan in each direction to compute the local mean coverage
    let scan_length = max(5000, 10 * max_read_len) as u32;
    // We want to fetch records from the region which is scan_length
    // before the start and after the end
    let start = locus.start.saturating_sub(scan_length);
    let end = locus.end + scan_length;

    let skip_length = 2 * max_read_len as u32; // To get rid of boundary bias

    // Find the coverage before locus.start
    bam.fetch(tid,
              start.saturating_sub(skip_length),
              locus.start + skip_length)
        .expect("Error fetching BAM file.");
    let mut total_coverage = 0;
    // pileup over all covered sites
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let dist_from_boundary = min(pileup.pos() - start, end - pileup.pos());
        if dist_from_boundary <= skip_length {
            continue;
        }
        for alignment in pileup.alignments() {
            let record = alignment.record();
            // Skip poor mappings & secondary alignments
            if record.mapq() < 10 || record.is_secondary() {
                continue;
            }
            total_coverage += 1;
        }
    }

    // Find the coverage after locus.end
    bam.fetch(tid,
              locus.end.saturating_sub(skip_length),
              end + skip_length)
        .expect("Error fetching BAM file.");
    // pileup over all covered sites
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let dist_from_boundary = min(pileup.pos() - start, end - pileup.pos());
        if dist_from_boundary <= skip_length {
            continue;
        }
        for alignment in pileup.alignments() {
            let record = alignment.record();
            // Skip poor mappings & secondary alignments
            if record.mapq() < 10 || record.is_secondary() {
                continue;
            }
            total_coverage += 1;
        }
    }

    let mean_coverage = total_coverage / (locus.start - start + end - locus.end);

    if mean_coverage < min_mean_coverage {
        return *event;
    }

    // Find the merged event Range
    let mut start = std::usize::MAX;
    let mut end = 0;
    for range in &del_ranges {
        start = min(start, range.start);
        end = max(start, range.end);
    }
    let event_range = Range { start, end };

    // Find phased coverage in the event range
    bam.fetch(tid,
              event_range.start.saturating_sub(max_read_len) as u32,
              (event_range.end + max_read_len) as u32)
        .expect("Error fetching BAM file.");

    let mut phased_coverage = vec![vec![0,0,0]; (event_range.end-event_range.start) as usize];
    for p in bam.pileup() {
        let pileup = p.unwrap();
        if (pileup.pos() < event_range.start as u32) || (pileup.pos() >= event_range.end as u32) {
            continue;
        }
        let idx = pileup.pos() as usize - event_range.start;
        for alignment in pileup.alignments() {
            let record = alignment.record();
            // Skip poor mappings & secondary alignments
            if record.mapq() < 10 || record.is_secondary() {
                continue;
            }
            let hap = call::get_haplotype(&record).unwrap_or(0);
            phased_coverage[idx][hap as usize] += 1;
        }
    }

    let mut hyp_scores = Vec::new();

    for &(h1, h2) in &hypotheses {
        let mut hyp_score = 0.0f64;
        for (i, c) in phased_coverage.iter().enumerate() {
            let mut pos_score = 0.0f64;
            let base_pos = event_range.start + i;
            for k in 0..(c[0] + 1) {
                let p1 = base_probability(base_pos, &del_ranges[h1], c[1] + k, mean_coverage);
                let p2 =
                    base_probability(base_pos, &del_ranges[h2], c[2] + c[0] - k, mean_coverage);
                pos_score += combinations(u64::from(c[0]), u64::from(k)) * p1 * p2
            }
            hyp_score += pos_score.ln() - f64::from(c[0]) * 2.0f64.ln();
        }
        if hyp_score.is_nan() || hyp_score.is_infinite() {
            // combinations can get big
            return *event;
        }
        hyp_scores.push((hyp_score, (h1, h2)));
    }

    hyp_scores.sort_by(|x, y| x.partial_cmp(y).unwrap());
    info!("coverage validation scores: {:?}", hyp_scores);

    let (h0, h1) = hyp_scores[hyp_scores.len() - 1].1;
    let mut haps = vec![event.0, event.1];
    haps.retain(|&x| x.is_some());
    let e0 = if h0 == 0 { None } else { haps[h0 - 1] };
    let e1 = if h1 == 0 { None } else { haps[h1 - 1] };
    (e0, e1)

}

// Model for P(Base Deletion | Coverage)
// Exponential with probability of 0.99 at 0 coverage and 0.02 at mean coverage
fn p_base_del_coverage(coverage: u32, mean_coverage: u32) -> f64 {
    let x = f64::from(coverage) / f64::from(mean_coverage);
    0.99f64 * (-3.9f64 * x).exp()
}

// Given a base position, deletion range and coverage, returns
// P(Base Deletion | Coverage) if the base position is withing deletion range
// Otherwise returns P(No Base Deletion | Coverage) = 1 - P(Base Deletion | Coverage)
fn base_probability(base_pos: usize,
                    del_range: &Range<usize>,
                    coverage: u32,
                    mean_coverage: u32)
                    -> f64 {
    if (base_pos >= del_range.start) && (base_pos < del_range.end) {
        p_base_del_coverage(coverage, mean_coverage)
    } else {
        1.0f64 - p_base_del_coverage(coverage, mean_coverage)
    }
}

fn score_hap(read_pairs: &Vec<HapReadPair>, hap: &[u8], aligner: &mut Aligner) -> Vec<i32> {
    let hap_rc = alphabets::dna::revcomp(hap);

    let mut scores = vec![];
    for rp in read_pairs.iter() {
        let sc = score_hap_read_pair(rp, hap, &hap_rc, aligner);
        scores.push(sc);
    }

    scores
}

pub fn validate(bam: &mut IndexedReader,
                locus: &Locus,
                _ref_a: &[u8],
                ref_start: usize,
                event: (Option<((usize, usize), usize)>, Option<((usize, usize), usize)>))
                -> (Option<((usize, usize), usize)>, Option<((usize, usize), usize)>) {


    let aln_call = validate_using_alignment(bam, locus, _ref_a, ref_start, &event);
    if aln_call == (None, None) {
        return aln_call;
    }

    validate_using_coverage(bam, locus, &aln_call)

}
pub fn validate_using_alignment
(bam: &mut IndexedReader,
 locus: &Locus,
 _ref_a: &[u8],
 ref_start: usize,
 event: &(Option<((usize, usize), usize)>, Option<((usize, usize), usize)>))
 -> (Option<((usize, usize), usize)>, Option<((usize, usize), usize)>) {

    use bio::alignment::pairwise::{Scoring, banded};

    // Need to make sure that the sparse DP gets the right match score
    let score = Scoring::from_scores(-6, -1, 1, -4).xclip(-5).yclip(0);
    let mut aligner = banded::Aligner::with_scoring(score, 12, 6);

    let ref_a = Vec::from(_ref_a);
    let rps = get_read_pairs(bam, locus);

    // Appy the edit e to reference r
    let apply = |e: ((usize, usize), usize)| {
        let ((ev_len, ev_pos), _) = e;
        let ev_start = ev_pos - ref_start;
        let mut edit = vec![];

        edit.extend(&_ref_a[0..ev_start]);
        edit.extend(&_ref_a[(ev_start + ev_len).._ref_a.len()]);
        edit
    };

    let (haplotypes, hypotheses) = match *event {
        (None, None) => panic!("don't validate empty events"),
        (None, Some(e)) => (vec![ref_a, apply(e)], vec![(0, 0), (0, 1), (1, 1)]),
        (Some(e), None) => (vec![ref_a, apply(e)], vec![(0, 0), (1, 0), (1, 1)]),
        (Some(e1), Some(e2)) if e1 == e2 => {
            (vec![ref_a, apply(e1)], vec![(0, 0), (0, 1), (1, 0), (1, 1)])
        }
        (Some(e1), Some(e2)) => {
            (vec![ref_a, apply(e1), apply(e2)], vec![(0, 0), (1, 0), (0, 2), (1, 2)])
        }
    };

    let mut hap_read_scores = vec![];
    for h in haplotypes {
        hap_read_scores.push(score_hap(&rps, &h, &mut aligner));
    }

    /*
    println!("hyps: {:?}", hypotheses);
    for i in 0..rps.len() {
        print!("Hap {} {} -> \t", rps[i].hap.unwrap_or(0), rps[i].seq1.len());
        for h in 0..hap_read_scores.len() {
            print!("{:?}\t", hap_read_scores[h][i]);
        }
        println!("");
    }
    */

    let mut hyp_scores = vec![];

    for &(h1, h2) in &hypotheses {
        // Score a hypothesis (a ordered pair of haplotypes), against the read-pair scores

        let mut hyp_score = 0.0;

        for (idx, rp) in rps.iter().enumerate() {

            let p = if let Some(q) = rp.hap_qual {
                1.0f64 - 10.0f64.powf(-f64::from(q) / 10.0f64)
            } else {
                0.999_99f64
            };

            let (ln_p_h1, ln_p_h2) = match rp.hap {
                None => (0.5f64.ln(), 0.5f64.ln()),
                Some(1) => (p.ln(), (1.0f64 - p).ln()),
                Some(2) => ((1.0f64 - p).ln(), p.ln()),
                _ => panic!("bad hap value"),
            };

            let _h1 = ln_p_h1 + f64::from(hap_read_scores[h1][idx]);
            let _h2 = ln_p_h2 + f64::from(hap_read_scores[h2][idx]);

            let read_len = rp.seq1.len() +
                match rp.seq2 {
                    Some(ref seq) => seq.len(),
                    None => 0,
                };
            let noise = 0.3f64 * (read_len as f64);

            let read_score = LogProb::ln_sum_exp(&[LogProb(_h1), LogProb(_h2), LogProb(noise)]).0;
            hyp_score += read_score;
        }

        hyp_scores.push(hyp_score);
    }

    let mut hyp_scores: Vec<_> = hyp_scores.iter().zip(hypotheses).collect();
    hyp_scores.sort_by(|x, y| x.partial_cmp(y).unwrap());

    info!("alignment validation scores: {:?}", hyp_scores);
    let (h0, h1) = hyp_scores[hyp_scores.len() - 1].1;
    let mut haps = vec![event.0, event.1];
    haps.retain(|&x| x.is_some());
    let e0 = if h0 == 0 { None } else { haps[h0 - 1] };
    let e1 = if h1 == 0 { None } else { haps[h1 - 1] };
    (e0, e1)
}


fn write_hits<P: AsRef<Path>>(p: P, read: &[u8], target: &[u8]) {
    //println!("read: {}, target:{}", read.len(), target.len());
    let hits = sparse::find_kmer_matches(read, target, 8);
    write_table(&hits, p);
}

fn score_seq_pacbio(read: &[u8], target: &[u8]) -> Alignment {
    let score = |a: u8, b: u8| if a == b { 2i32 } else { -7i32 };

    let mut banded_aligner =
        banded::Aligner::with_capacity(read.len(), target.len(), -100, -5, &score, 8, 50);
    banded_aligner.local(read, target)

}



fn score_del(del: (String, u32, u32),
             bam: &mut IndexedReader,
             fa: &mut fasta::IndexedReader<File>,
             debug: bool,
             from_asm: bool)
             -> Vec<SvObs> {

    let sz = del.2 - del.1;

    let pad = if from_asm { 1000 + sz } else { 300 + sz / 3 };

    let locus = Locus {
        chrom: del.0,
        start: max(0, del.1.saturating_sub(pad)),
        end: del.2 + pad,
    };

    if sz > 20_000 {
        return Vec::new();
    }

    if debug {
        println!("Locus: {:?}", locus);
    }

    let (ref_seq, _) = call::read_locus(fa, &locus, 0, 0);

    let mut alt_seq = Vec::new();
    alt_seq.extend_from_slice(&ref_seq[..(del.1 - locus.start) as usize]);
    alt_seq.extend_from_slice(&ref_seq[(del.2 - locus.start) as usize..]);

    if debug {
        println!("ref: {}", String::from_utf8_lossy(&ref_seq));
        println!("alt: {}", String::from_utf8_lossy(&alt_seq));
    }


    let mut scores = Vec::new();

    let tid = bam.header.tid(locus.chrom.as_bytes()).unwrap();
    bam.fetch(tid, locus.start, locus.end)
        .expect("Error fetching BAM file.");

    for _rec in bam.records() {
        let rec = _rec.ok().unwrap();
        let r = rec.seq().as_bytes();

        if rec.mapq() < 10 || r.len() < 50 {
            continue;
        }

        let read = if from_asm {
            call::get_seq_bounded(&rec, &locus)
        } else {
            // For checking SN contigs against events, allow secondary / supplementary.
            if rec.is_secondary() || rec.is_supplementary() {
                continue;
            }
            rec.seq().as_bytes()
        };

        if read.len() > 30_000 {
            continue;
        }
        let qname = str::from_utf8(rec.qname()).unwrap();

        let ref_aln = score_seq_pacbio(&read, &ref_seq);
        let alt_aln = score_seq_pacbio(&read, &alt_seq);

        if debug {
            write_hits(format!("hits_{}_ref.tsv", qname), &read, &ref_seq);
            write_hits(format!("hits_{}_alt.tsv", qname), &read, &alt_seq);

            //write_table(&ref_aln.path(), format!("aln_{}_ref.tsv", qname));
            //write_table(&alt_aln.path(), format!("aln_{}_alt.tsv", qname));
        }

        //let diff_score = (alt_aln.score as f32 - ref_aln.score as f32) / (del.2 as f32);

        let obs = SvObs {
            chrom: locus.chrom.clone(),
            start: locus.start,
            end: locus.end,
            is_reverse: rec.is_reverse(),
            qname: qname.to_string(),
            read_length: read.len() as u32,
            ref_score: ref_aln.score,
            alt_score: alt_aln.score,
        };

        if debug {
            println!("{:?}", obs);
        }
        scores.push(obs);
    }

    scores
}


pub fn write_table<P: AsRef<Path>, T: Serialize>(data: &Vec<T>, path: P) {
    let mut w = csv::WriterBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    for d in data {
        w.serialize(d).unwrap();
    }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
struct SvObs {
    chrom: String,
    start: u32,
    end: u32,
    is_reverse: bool,
    qname: String,
    read_length: u32,
    ref_score: i32,
    alt_score: i32,
}

#[derive(Serialize, Deserialize)]
pub struct ValObs {
    name: String,
    total_reads: usize,
    alt_support_reads: usize,
    alt_support_score: i32,
    ref_support_reads: usize,
    ref_support_score: i32,
    uninformative_reads: usize,
    uninformative_score: i32,
}



pub fn go_one(args: &Args) {
    let mut fa = fasta::IndexedReader::from_file(&args.arg_fasta).expect("error opening fa");
    let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
        .expect("Error opening BAM file");

    let locus: Locus = FromStr::from_str(&args.arg_locus.clone().unwrap()).unwrap();

    let r = score_del((locus.chrom, locus.start, locus.end),
                      &mut bam,
                      &mut fa,
                      args.flag_debug,
                      args.flag_asm);
    write_table(&r, "validation_reads.tsv");
    if args.flag_asm {
        for ev in &r {
            println!("{:?}", ev);
        }
    }

    let mut nref = 0;
    let mut nalt = 0;
    let mut sref = 0;
    let mut salt = 0;
    for o in r {
        if o.ref_score > o.alt_score {
            nref += 1;
            sref += o.ref_score - o.alt_score;
        }
        if o.alt_score > o.ref_score {
            nalt += 1;
            salt += o.alt_score - o.ref_score;
        }
    }
    println!("N ref: {}, Support ref:{}, N alt:{}, Support alt:{}",
             nref,
             sref,
             nalt,
             salt);
}

pub fn go_bedpe(args: &Args) {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&args.arg_bed)
        .unwrap();
    let bedpe_rows = rdr.deserialize()
        .collect::<csv::Result<Vec<BedpeRow>>>()
        .unwrap();
    let sv_validations = process_chunk(&bedpe_rows, args);

    write_table(&sv_validations, "sv_validations.tsv");
    println!("done validation");
}


pub fn process_chunk(bedpe_rows: &[BedpeRow], args: &Args) -> Vec<ValObs> {

    let mut fa = fasta::IndexedReader::from_file(&args.arg_fasta).expect("error opening fa");
    let mut bam = bam::IndexedReader::from_path(&args.arg_bam)
        .expect("Error opening BAM file");

    let mut targets = HashSet::new();
    for n in bam.header.target_names() {
        targets.insert(String::from_utf8(n.to_vec()).unwrap());
    }

    // load bedpe

    let mut sv_validations = Vec::new();

    for bedpe_entry in bedpe_rows {
        println!("{:?}", bedpe_entry);

        let start = (bedpe_entry.start1 + bedpe_entry.end1) / 2;
        let end = (bedpe_entry.start2 + bedpe_entry.end2) / 2;
        let event_sz = (end - start) as i32;

        if !targets.contains(&bedpe_entry.chrom1) {
            continue;
        }

        let r = score_del((bedpe_entry.chrom1.clone(), start as u32, end as u32),
                          &mut bam,
                          &mut fa,
                          args.flag_debug,
                          args.flag_asm);

        let mut sv_validation = ValObs {
            name: bedpe_entry.name.clone(),
            total_reads: r.len(),
            alt_support_reads: 0,
            alt_support_score: 0,
            ref_support_reads: 0,
            ref_support_score: 0,
            uninformative_reads: 0,
            uninformative_score: 0,
        };

        for obs in r {
            let diff = obs.alt_score - obs.ref_score;

            if diff * 10 > obs.ref_score || diff * 3 > event_sz {
                sv_validation.alt_support_reads += 1;
                sv_validation.alt_support_score += diff;
            } else if -diff * 10 > obs.alt_score || -diff * 3 > event_sz {
                sv_validation.ref_support_reads += 1;
                sv_validation.ref_support_score += diff;
            } else {
                sv_validation.uninformative_reads += 1;
                sv_validation.uninformative_score += diff;
            }
        }

        sv_validations.push(sv_validation);
    }

    sv_validations
}