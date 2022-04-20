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
        sample_idx: usize,
    ) -> (usize, f64);
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
        sample_idx: usize,
    ) -> (usize, f64) {
        match gene.strand() {
            Some(strand) => {

                let contig_name = gene.seqname();
                let rid = if let Some(rid) = VariantContext::get_contig_vcf_tid(variants.header(), contig_name.as_bytes()) {
                    rid
                } else {
                    // no variants on this contig so skip
                    return (0, 1.0)
                };

                reference_reader.fetch_contig_from_reference_by_tid(rid as usize, ref_idx);

                // bio::gff documentation says start and end positions are 1-based, so we minus 1
                // Additionally, end position is non-inclusive so do minus 1
                let start = *gene.start() as usize - 1;
                let end = *gene.end() as usize - 1;

                // fetch variants in this window
                variants.fetch(rid, start as u64, Some(end as u64));

                let frame: usize = match gene.frame().parse() {
                    Ok(frame_val) => frame_val,
                    Err(_) => 0,
                };
                let gene_sequence = &reference_reader.current_sequence[start..=end];
                debug!("Gene Seq {:?}", String::from_utf8_lossy(gene_sequence));
                let codon_sequence = get_codons(&gene_sequence, frame, strand);
                debug!("Codon Sequence {:?}", codon_sequence);

                // Calculate N and S
                let mut big_n: f64 = 0.0;
                let mut big_s: f64 = 0.0;
                for codon in codon_sequence.iter() {
                    if std::str::from_utf8(codon)
                        .expect("Unable to interpret codon")
                        .contains('N')
                        || codon.len() != 3
                    {
                        continue;
                    } else {
                        let n = self.ns_sites.get(codon).unwrap();
                        big_n += n;
                        big_s += 3.0 - n;
                    }
                }

                debug!("getting ns_sites N {} S {}", big_n, big_s);

                // Create Nd and Sd values
                let mut big_nd: f64 = 0.0;
                let mut big_sd: f64 = 0.0;

                // dN/dS calculations when using NGS reads outlined here:
                // http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
                // Note, we don't normalize for depth here and instead just use Jukes-Cantor model
                let mut codon: Vec<u8> = vec![];
                let mut new_codons: Vec<Vec<u8>> = vec![];
                let mut positionals = 0;
                let mut total_variants = 0;
                let mut gene_cursor = 0; // index inside gene, i.e, position of variant
                let mut reference_cursor = start; // index of reference sequence
                let mut frameshifts = 0;
                for record in variants.records().into_iter() {
                    match record {
                        Ok(mut record) => {
                            match VariantContext::from_vcf_record(&mut record, true) {
                                Some(context) => {
                                    // gained bases is the difference between current and previous
                                    // position
                                    let gained_bases = context.loc.start - gene_cursor;
                                    // update previous position to current position
                                    gene_cursor = context.loc.start;
                                    // add gained bases to reference cursor
                                    reference_cursor += gained_bases;

                                    // index of current codon in gene
                                    // number of codons is gene size divided by three
                                    let codon_idx = gene_cursor / 3 as usize;
                                    // index of current position in the current codon
                                    // all codons have size == 3
                                    let codon_cursor = gene_cursor % 3;

                                    if std::str::from_utf8(codon.as_slice())
                                        .expect("Unable to interpret codon")
                                        .contains('N')
                                    {
                                        continue;
                                    }

                                    if codon_cursor == 0 || new_codons.len() == 0 {
                                        for new_codon in new_codons.iter_mut() {
                                            if (codon.len() == 3) && (new_codon.len() == 3) && (codon != *new_codon)
                                            {
                                                // get indices of different locations
                                                let mut pos = 0 as usize;
                                                let mut diffs = vec![];
                                                for (c1, c2) in codon.iter().zip(new_codon.iter()) {
                                                    if c1 != c2 {
                                                        diffs.push(pos);
                                                    }
                                                    pos += 1;
                                                }
                                                total_variants += diffs.len();
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
                                                positionals += permutations.len();
                                                for permutation in permutations.iter() {
                                                    let mut shifting = codon.clone();
                                                    let mut old_shift;
                                                    for pos in permutation {
                                                        // Check if one amino acid change causes an syn or non-syn
                                                        old_shift = shifting.clone();
                                                        shifting[**pos] = new_codon[**pos];
                                                        debug!("Old shift {:?}, new {:?}", old_shift, shifting);
                                                        if self.aminos[&old_shift] != self.aminos[&shifting] {
                                                            ns += 1;
                                                        } else {
                                                            ss += 1;
                                                        }
                                                    }
                                                }
                                                let nd = ns as f64 / permutations.len() as f64;
                                                let sd = ss as f64 / permutations.len() as f64;
                                                big_nd += nd;
                                                big_sd += sd;
                                            }
                                        }
                                        // begin working on new codon
                                        debug!(
                                            "Codon idx {} codonds {} gene length {} new_codons {:?}",
                                            codon_idx,
                                            codon_sequence.len(),
                                            gene_sequence.len(),
                                            new_codons
                                        );
                                        codon = codon_sequence[codon_idx].clone();
                                        if codon.len() != 3 {
                                            continue;
                                        }
                                        new_codons = vec![codon.clone()];
                                    }

                                    let ref_allele = context.get_reference();
                                    let mut snp_count = 0;

                                    // iterate through non reference alleles
                                    // if those alleles are present in this sample then
                                    // increment appropriate values
                                    for (allele_index, allele) in context.get_alternate_alleles_with_index() {
                                        if allele.bases.len() > 1 || allele.bases.len() != ref_allele.bases.len() {
                                            if context.get_genotypes().genotypes()[sample_idx].ad[allele_index] > 0 {
                                                frameshifts += 1;
                                            }
                                            continue
                                        }

                                        if context.get_genotypes().genotypes()[sample_idx].ad[allele_index] > 0 {
                                            snp_count += 1;
                                            if snp_count > 1 {
                                                // Create a copy of codon up to this point
                                                // Not sure if reusing previous variants is bad, but
                                                // not doing so can cause randomness to dN/dS values

                                                new_codons.push(codon.clone());

                                                new_codons[snp_count][codon_cursor] = allele.bases[0];

                                                debug!("multi snp codon {:?}", new_codons);
                                            } else {
                                                for var_idx in 0..new_codons.len() {
                                                    new_codons[var_idx][codon_cursor] = allele.bases[0];
                                                }
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

                debug!(
                    "Nd {} N {}, Sd {} S {} total permutations {} variants {}",
                    big_nd, big_n, big_sd, big_s, positionals, total_variants
                );
                let mut pn = big_nd / big_n;
                let mut ps = big_sd / big_s;
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
                if dnds < 0.0 {
                    dnds = 0.0
                }

                if dnds.is_nan() {
                    dnds = 1.
                }
                return (frameshifts, dnds);
            }
            _ => return (0, 1.0),
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

#[allow(unused)]
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

    #[test]
    fn test_dnds() {
        let mut codon_table = CodonTable::setup();
        codon_table.get_codon_table(11);

        let ref_sequence = "ATGAAACCCGGGTTTTAA".as_bytes().to_vec();
        let sample_count = 2;

        let mut gene_records = gff::Reader::from_file("tests/data/dnds.gff", gff::GffType::GFF3)
            .expect("Incorrect file path");
        let mut variant_abundances: HashMap<i64, HashMap<Variant, Base>> = HashMap::new();
        variant_abundances.insert(13, HashMap::new());
        variant_abundances.insert(14, HashMap::new());
        let hash = variant_abundances.entry(7).or_insert(HashMap::new());

        // Create fake variants
        let var_1 = create_base(&ref_sequence, "G".bytes().nth(0).unwrap(), 7, 2);
        let var_2 = create_base(&ref_sequence, "C".bytes().nth(0).unwrap(), 11, 2);
        let var_3 = create_base(&ref_sequence, "A".bytes().nth(0).unwrap(), 13, 2);
        let var_4 = create_base(&ref_sequence, "C".bytes().nth(0).unwrap(), 14, 2);

        hash.insert(var_1.variant.clone(), var_1);

        let hash = variant_abundances.entry(11).or_insert(HashMap::new());
        hash.insert(var_2.variant.clone(), var_2);

        let hash = variant_abundances.entry(13).or_insert(HashMap::new());
        hash.insert(var_3.variant.clone(), var_3);

        let hash = variant_abundances.entry(14).or_insert(HashMap::new());
        hash.insert(var_4.variant.clone(), var_4);

        for gene_record in gene_records.records() {
            let gene_record = gene_record.unwrap();

            let dnds =
                codon_table.find_mutations(&gene_record, &variant_abundances, &ref_sequence, 0);
            assert_eq!(format!("{:.4}", dnds), format!("{}", 0.1247));
        }
    }
}
