use bio::alphabets::dna;
use bio_types::strand::Strand;
use itertools::{izip, Itertools};
use model::variants::{Base, Variant};
use std::collections::HashMap;
use utils::utils::{mean, std_deviation};
use model::variant_context::VariantContext;
use model::variant_context_utils::VariantContextUtils;
use rust_htslib::bcf::Read;
use reference::reference_reader::ReferenceReader;
use rand::seq::index::sample;

#[allow(dead_code)]
pub struct GeneInfo {
    name: String,
    locations: HashMap<u32, (u64, u64)>,
    coverages: HashMap<u32, f64>,
}

pub struct CodonTable {
    aminos: HashMap<Vec<u8>, char>,
    starts: HashMap<Vec<u8>, char>,
    ns_sites: HashMap<Vec<u8>, f64>,
}

pub struct NCBITable {
    aas: String,
    starts: String,
    base1: String,
    base2: String,
    base3: String,
}

impl NCBITable {
    // get translation tables in NCBI format
    // Kind of lazy storing and then converting every time but would take way too much time
    // to write out each table into CodonTable format by hand
    fn get_translation_table(table_id: usize) -> NCBITable {
        match table_id {
            1 => NCBITable {
                aas: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".to_owned(),
                starts: "---M------**--*----M---------------M----------------------------"
                    .to_owned(),
                base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
                    .to_owned(),
                base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
                    .to_owned(),
                base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
                    .to_owned(),
            },
            11 => NCBITable {
                aas: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".to_owned(),
                starts: "---M------**--*----M------------MMMM---------------M------------"
                    .to_owned(),
                base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
                    .to_owned(),
                base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
                    .to_owned(),
                base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
                    .to_owned(),
            },
            _ => {
                panic!("Translation table {} not yet implemented", table_id);
            }
        }
    }
}

impl CodonTable {
    pub fn setup() -> CodonTable {
        CodonTable {
            aminos: HashMap::new(),
            starts: HashMap::new(),
            ns_sites: HashMap::new(),
        }
    }
}

pub trait Translations {
    fn get_codon_table(&mut self, table_id: usize);
    fn find_mutations(
        &self,
        gene: &bio::io::gff::Record,
        variants: &mut rust_htslib::bcf::IndexedReader,
        reference_reader: &mut ReferenceReader,
        ref_idx: usize,
        n_samples: usize,
        qual_by_depth_filter: f64,
        qual_threshold: f64,
        depth_per_sample_filter: i64
    ) -> (Vec<usize>, Vec<usize>, Vec<f64>);
    fn calculate_gene_coverage(
        &self,
        gene: &bio::io::gff::Record,
        depth_of_contig: &Vec<i32>,
    ) -> (f32, f32);
}

impl Translations for CodonTable {
    fn get_codon_table(&mut self, table_id: usize) {
        // Convert NCBI format to something more usable
        let ncbi_format = NCBITable::get_translation_table(table_id);
        let mut amino_hash = HashMap::new();
        let mut start_hash = HashMap::new();
        for (aa, s, b1, b2, b3) in izip!(
            ncbi_format.aas.as_bytes(),
            ncbi_format.starts.as_bytes(),
            ncbi_format.base1.as_bytes(),
            ncbi_format.base2.as_bytes(),
            ncbi_format.base3.as_bytes()
        ) {
            let codon = vec![*b1, *b2, *b3];
            amino_hash.insert(codon.clone(), *aa as char);
            start_hash.insert(codon, *s as char);
        }
        debug!("Amino hash {:?}", amino_hash);
        self.aminos = amino_hash;
        self.starts = start_hash;

        let nucleotides: Vec<u8> = vec![65, 84, 67, 71];
        debug!("nucleotides: {} ", String::from_utf8_lossy(&nucleotides));
        for (codon, _aa) in self.aminos.iter() {
            let mut n = 0.0;
            for (pos, cod) in codon.iter().enumerate() {
                for nuc in nucleotides.iter() {
                    if cod == nuc {
                        // Ignore this nucleotide
                        continue;
                    } else {
                        let mut codon_shift = codon.clone();
                        // Change one position of the codon
                        codon_shift[pos] = nuc.clone();

                        if self.aminos[codon] != self.aminos[&codon_shift] {
                            // This change can cause a non-synonymous mutation
                            n += 1.0 / 3.0;
                        }
                    }
                }
            }
            self.ns_sites.insert(codon.clone(), n as f64);
        }
    }

    /// Finds all associate mutations within a gene region in the form of a gff record
    /// If there are associated variants in this gene attempts to calculate dN/dS ratios for
    /// the given sample
    /// Returns a tuple of the number of frameshifts and dN/dS ratio
    /// TODO: Refactor so calculates for all samples at once without having to re-read the variant
    ///       region each time.
    fn find_mutations(
        &self,
        gene: &bio::io::gff::Record,
        variants: &mut rust_htslib::bcf::IndexedReader,
        reference_reader: &mut ReferenceReader,
        ref_idx: usize,
        n_samples: usize,
        qual_by_depth_filter: f64,
        qual_threshold: f64,
        depth_per_sample_filter: i64
    ) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        match gene.strand() {
            Some(strand) => {

                let contig_name = format!("{}~{}", reference_reader.retrieve_reference_stem(ref_idx), gene.seqname()); // create concatenated contig name format
                let rid = if let Some(rid) = VariantContext::get_contig_vcf_tid(variants.header(), contig_name.as_bytes()) {
                    rid
                } else {
                    // no variants on this contig so skip
                    return (vec![0; n_samples], vec![0; n_samples], vec![1.0; n_samples])
                };

                reference_reader.fetch_contig_from_reference_by_contig_name(contig_name.as_bytes(), ref_idx);
                reference_reader.read_sequence_to_vec();
                // bio::gff documentation says start and end positions are 1-based, so we minus 1
                // Additionally, end position is non-inclusive so do minus 1
                let start = *gene.start() as usize - 1;
                let end = *gene.end() as usize - 1;
                debug!("Start {} End {}", start, end);
                // fetch variants in this window
                variants.fetch(rid, start as u64, Some(end as u64));

                // VariantContext::process_vcf_in_region()

                let frame: usize = match gene.frame().parse() {
                    Ok(frame_val) => frame_val,
                    Err(_) => 0,
                };
                let gene_sequence = &reference_reader.current_sequence[start..=end];
                debug!("Gene Seq {:?}", String::from_utf8_lossy(gene_sequence));
                let codon_sequence = get_codons(&gene_sequence, frame, strand);
                debug!("Codon Sequence {:?}", codon_sequence);

                // Calculate N and S
                let mut big_n: Vec<f64> = vec![0.0; n_samples];
                let mut big_s: Vec<f64> = vec![0.0; n_samples];
                for codon in codon_sequence.iter() {
                    if std::str::from_utf8(codon)
                        .expect("Unable to interpret codon")
                        .contains('N')
                        || codon.len() != 3
                    {
                        continue;
                    } else {
                        let n = self.ns_sites.get(codon).unwrap();
                        for sample_idx in 0..n_samples {
                            big_n[sample_idx] += n;
                            big_s[sample_idx] += 3.0 - n;
                        }
                    }
                }

                debug!("getting ns_sites N {:?} S {:?}", &big_n, &big_s);

                // Create Nd and Sd values
                let mut big_nd: Vec<f64> = vec![0.0; n_samples];
                let mut big_sd: Vec<f64> = vec![0.0; n_samples];

                // dN/dS calculations when using NGS reads outlined here:
                // http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
                // Note, we don't normalize for depth here and instead just use Jukes-Cantor model
                let mut codon: &Vec<u8> = &codon_sequence[0];
                let mut previous_codon: &Vec<u8> = &codon_sequence[0];
                let mut new_codons: Vec<Vec<Vec<u8>>> = vec![vec![]; n_samples];
                let mut positionals = vec![0; n_samples];
                let mut total_variants = vec![0; n_samples];
                let mut gene_cursor = 0; // index inside gene, i.e, position of variant
                let mut reference_cursor = start; // index of reference sequence
                let mut frameshifts = vec![0; n_samples];
                let mut snps = vec![0; n_samples];
                let mut old_codon_idx = None;
                for record in variants.records().into_iter() {
                    match record {
                        Ok(mut record) => {
                            match VariantContext::from_vcf_record(&mut record, true) {
                                Some(mut context) => {
                                    let passes = VariantContextUtils::passes_thresholds(
                                        &mut context, qual_by_depth_filter, qual_threshold
                                    );

                                    if !passes {
                                        continue
                                    }

                                    // gained bases is the difference between current and previous
                                    // position
                                    let gained_bases = context.loc.start.saturating_sub(reference_cursor);
                                    // update previous position to current position
                                    gene_cursor += gained_bases;
                                    // add gained bases to reference cursor
                                    reference_cursor = context.loc.start;

                                    // index of current codon in gene
                                    // number of codons is gene size divided by three
                                    let codon_idx = gene_cursor / 3 as usize;
                                    let process_previous_codon = match old_codon_idx {
                                        Some(old_idx) => {
                                            if old_idx != codon_idx {
                                                old_codon_idx = Some(codon_idx);
                                                previous_codon = codon;
                                                true
                                            } else {
                                                false
                                            }
                                        },
                                        None => {
                                            old_codon_idx = Some(codon_idx);
                                            false
                                        }
                                    };


                                    // index of current position in the current codon
                                    // all codons have size == 3
                                    let codon_cursor = gene_cursor % 3;
                                    debug!(
                                        "reference cursor {} gained_bases {} Gene cursor {} codon idx {} codon cursor {}",
                                        reference_cursor, gained_bases, gene_cursor, codon_idx, codon_cursor
                                    );

                                    if codon_idx >= codon_sequence.len() {
                                        continue
                                    }
                                    codon = &codon_sequence[codon_idx];


                                    if std::str::from_utf8(codon.as_slice())
                                        .expect("Unable to interpret codon")
                                        .contains('N') || codon.len() != 3
                                    {
                                        continue;
                                    }

                                    for sample_idx in 0..n_samples {
                                        if process_previous_codon {
                                            for new_codon in new_codons[sample_idx].iter() {
                                                debug!("Positionals start {:?}", &positionals);
                                                debug!("Sample {} new_codon {:?} ref_codon {:?}", sample_idx, new_codon, &previous_codon);
                                                if (previous_codon.len() == 3) && (new_codon.len() == 3) && (previous_codon != new_codon)
                                                {
                                                    // get indices of different locations
                                                    let mut pos = 0 as usize;
                                                    let mut diffs = vec![];
                                                    for (c1, c2) in previous_codon.iter().zip(new_codon.iter()) {
                                                        if c1 != c2 {
                                                            diffs.push(pos);
                                                        }
                                                        pos += 1;
                                                    }
                                                    total_variants[sample_idx] += diffs.len();
                                                    // get permuations of positions
                                                    let permutations: Vec<Vec<&usize>> =
                                                        diffs.iter().permutations(diffs.len()).collect();

                                                    // calculate synonymous and non-synonymous for each permutation
                                                    let mut ns = 0;
                                                    let mut ss = 0;
                                                    debug!(
                                                        "positional difference {:?} permutations {:?}",
                                                        diffs,
                                                        permutations.len()
                                                    );
                                                    positionals[sample_idx] += permutations.len();
                                                    for permutation in permutations.iter() {
                                                        let mut shifting = codon.clone();
                                                        let mut old_shift;
                                                        for pos in permutation {
                                                            // Check if one amino acid change causes an syn or non-syn
                                                            old_shift = shifting.clone();
                                                            shifting[**pos] = new_codon[**pos];
                                                            // debug!("Old shift {:?}, new {:?}", old_shift, shifting);
                                                            if self.aminos[&old_shift] != self.aminos[&shifting] {
                                                                ns += 1;
                                                            } else {
                                                                ss += 1;
                                                            }
                                                        }
                                                    }
                                                    let nd = ns as f64 / permutations.len() as f64;
                                                    let sd = ss as f64 / permutations.len() as f64;
                                                    big_nd[sample_idx] += nd;
                                                    big_sd[sample_idx] += sd;
                                                }
                                                debug!("Positionals end {:?}", &positionals);
                                                debug!(
                                                    "Codon idx {} gene cursor {} codons {} gene length {} new_codons {:?}",
                                                    codon_idx,
                                                    gene_cursor,
                                                    codon_sequence.len(),
                                                    gene_sequence.len(),
                                                    new_codons[sample_idx]
                                                );
                                            }
                                            // begin working on new codon
                                            new_codons[sample_idx] = vec![codon.clone()];
                                        }

                                        let mut which_are_present =
                                            context.alleles_present_in_sample(
                                                sample_idx, depth_per_sample_filter);

                                        if !which_are_present[1..].iter().any(|v| *v) {
                                            continue // no alt alleles are present
                                        }

                                        let ref_allele = context.get_reference();
                                        let mut snp_count = 0;

                                        // iterate through non reference alleles
                                        // if those alleles are present in this sample then
                                        // increment appropriate values
                                        for (allele_index, allele) in context.get_alternate_alleles_with_index() {
                                            if allele.bases.len() > 1 || allele.bases.len() != ref_allele.bases.len() {
                                                if which_are_present[allele_index] {
                                                    frameshifts[sample_idx] += 1;
                                                }
                                                continue
                                            }

                                            if which_are_present[allele_index] {
                                                snps[sample_idx] += 1;
                                                if snp_count >= 1 {
                                                    // Create a copy of codon up to this point
                                                    // Not sure if reusing previous variants is bad, but
                                                    // not doing so can cause randomness to dN/dS values
                                                    debug!("before pushing {:?}", new_codons[sample_idx]);
                                                    new_codons[sample_idx].push(codon.clone());
                                                    debug!(
                                                        "Gene cursor {} sample idx {} snp count {} codon cursor {}",
                                                        gene_cursor, sample_idx, snp_count, codon_cursor
                                                    );
                                                    debug!("new_codons length {}", new_codons[sample_idx].len());
                                                    new_codons[sample_idx][snp_count][codon_cursor] = allele.bases[0];

                                                    debug!("multi snp codon {:?}", new_codons[sample_idx]);
                                                } else {
                                                    if new_codons[sample_idx].len() == 0 {
                                                        new_codons[sample_idx] = vec![codon.clone()];
                                                    }
                                                    for var_idx in 0..new_codons[sample_idx].len() {
                                                        debug!("s {} v {} c {}. Size {}", sample_idx, var_idx, codon_cursor, new_codons[sample_idx].len());
                                                        new_codons[sample_idx][var_idx][codon_cursor] = allele.bases[0];
                                                    }
                                                }
                                                snp_count += 1;
                                            }
                                        }
                                    }
                                },
                                None => {
                                    continue
                                }
                            }

                        },
                        Err(_) => {
                            // Skip error record
                            continue
                        }
                    }
                }

                let mut dnds_values = vec![1.0; n_samples];
                for sample_idx in 0..n_samples {
                    debug!(
                        "Nd {} N {}, Sd {} S {} total permutations {} variants {}",
                        big_nd[sample_idx], big_n[sample_idx], big_sd[sample_idx],
                        big_s[sample_idx], positionals[sample_idx], total_variants[sample_idx]
                    );
                    let mut pn = big_nd[sample_idx] / big_n[sample_idx];
                    let mut ps = big_sd[sample_idx] / big_s[sample_idx];
                    debug!("pn {} ps {}", pn, ps);
                    // Weirdly in the Jukes-Cantor model if pn or ps are 0.75 then the nat log does not resolve
                    // No one talks about this in the literature for some reason
                    if pn == 0.75 {
                        pn = 0.7499
                    }
                    if ps == 0.75 {
                        ps = 0.7499
                    }
                    let d_n = -(3.0 / 4.0) * (1.0 - (4.0 * pn) / 3.0).ln();
                    let d_s = -(3.0 / 4.0) * (1.0 - (4.0 * ps) / 3.0).ln();
                    debug!("dN {} dS {}", d_n, d_s);
                    let mut dnds = d_n / d_s;

                    // negative dnds values make no sense, but occur nonetheless
                    // Just make them 0.0
                    if dnds.is_sign_negative() {
                        dnds = 0.0
                    }

                    if dnds.is_nan() {
                        dnds = 1.
                    }
                    dnds_values[sample_idx] = dnds
                }

                return (snps, frameshifts, dnds_values);
            }
            _ => return (vec![0; n_samples], vec![0; n_samples], vec![1.0; n_samples]),
        }
    }

    fn calculate_gene_coverage(
        &self,
        gene: &bio::io::gff::Record,
        depth_of_contig: &Vec<i32>,
    ) -> (f32, f32) {
        let gene_start = *gene.start() as usize - 1;
        let gene_end = *gene.end() as usize - 1;
        let gene_depths = &depth_of_contig[gene_start..gene_end];
        let mean_cov = mean(gene_depths).unwrap_or(0.);
        let std = std_deviation(gene_depths).unwrap_or(0.);

        return (mean_cov, std);
    }
}

pub fn get_codons<'a>(sequence: &'a [u8], frame: usize, strandedness: Strand) -> Vec<Vec<u8>> {
    match strandedness {
        Strand::Forward | Strand::Unknown => sequence[0 + frame..]
            .chunks(3)
            .map(|chunk| chunk.to_vec())
            .collect::<Vec<Vec<u8>>>(),
        Strand::Reverse => {
            let rc = dna::revcomp(sequence);
            rc[0 + frame..]
                .chunks(3)
                .map(|chunk| chunk.to_vec())
                .collect::<Vec<Vec<u8>>>()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::gff;
    use bio::stats::LogProb;
    use model::variants;
    use std::collections::HashSet;

    fn create_base(ref_sequence: &Vec<u8>, var_char: u8, pos: i64, sample_count: usize) -> Base {
        Base {
            tid: 0,
            pos,
            refr: ref_sequence[pos as usize..(pos as usize + 1)].to_vec(),
            variant: Variant::SNV(var_char),
            depth: vec![5; sample_count],
            truedepth: vec![5; sample_count],
            totaldepth: vec![5; sample_count],
            genotypes: HashSet::new(),
            quals: vec![0.; sample_count],
            referencedepth: vec![0; sample_count],
            freq: vec![0.; sample_count],
            rel_abunds: vec![0.; sample_count],
            reads: HashSet::new(),
        }
    }

    #[test]
    fn test_floor_division() {
        assert_eq!(0 / 3 as usize, 0);
        assert_eq!(803 / 3 as usize, 267);
    }

}
