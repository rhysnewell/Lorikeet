use std::collections::{HashMap, HashSet};
use std::str;
use itertools::{izip, Itertools};
use bio::alphabets::dna;
use bio_types::strand;
use permutohedron::{Heap};


pub struct CodonTable {
    aminos: HashMap<Vec<u8>, char>,
    starts: HashMap<Vec<u8>, char>,
    ns_sites: HashMap<Vec<u8>, f32>,
}


pub struct NCBITable {
    AAs: String,
    Starts: String,
    Base1: String,
    Base2: String,
    Base3: String,
}


impl NCBITable {
    // get translation tables in NCBI format
    // Kind of lazy storing and then converting every time but would take way too much time
    // to write out each table into CodonTable format by hand
    fn get_translation_table(table_id: usize) -> NCBITable {
        match table_id {
            1 => {
                NCBITable {
                    AAs: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".to_owned(),
                    Starts: "---M------**--*----M---------------M----------------------------".to_owned(),
                    Base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".to_owned(),
                    Base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".to_owned(),
                    Base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".to_owned(),
                }
            },
            11 => {
                NCBITable {
                    AAs: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".to_owned(),
                    Starts: "---M------**--*----M------------MMMM---------------M------------".to_owned(),
                    Base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".to_owned(),
                    Base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".to_owned(),
                    Base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".to_owned(),
                }
            },
            _ => {
                panic!("Translation table {} not yet implemented", table_id);
            },
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
    fn find_mutations(&self,
                      gene: &bio::io::gff::Record,
                      variant_abundances: &Vec<HashMap<String, f32>>,
                      ref_sequence: &Vec<u8>,
                      depth: &Vec<usize>);
}

impl Translations for CodonTable {
    fn get_codon_table(&mut self, table_id: usize) {
        // Convert NCBI format to something more usable
        let ncbi_format = NCBITable::get_translation_table(table_id);
        let mut amino_hash = HashMap::new();
        let mut start_hash = HashMap::new();
        for (aa, s, b1, b2, b3) in izip!(ncbi_format.AAs.as_bytes(), ncbi_format.Starts.as_bytes(),
                                        ncbi_format.Base1.as_bytes(), ncbi_format.Base2.as_bytes(),
                                        ncbi_format.Base3.as_bytes()) {
            let mut codon = vec!(*b1, *b2, *b3);
            amino_hash.insert(codon.clone(), *aa as char);
            start_hash.insert(codon, *s as char);
        }
        debug!("Amino hash {:?}", amino_hash);
        self.aminos = amino_hash;
        self.starts = start_hash;

        let nucleotides: Vec<u8> = vec!("A" as u8, "T" as u8, "C" as u8, "G" as u8);
        for (codon, aa) in self.aminos{
            let n = 0.0;
            for (pos, cod) in codon.iter().enumerate() {
                for nuc in nucleotides {
                    if *cod == nuc {
                        // Ignore this nucleotide
                        continue
                    } else {
                        let mut codon_shift = codon.clone();
                        // Change one position of the codon
                        codon_shift[pos] = nuc;

                        if self.aminos[&codon] != self.aminos[&codon_shift] {
                            // This change can cause a non-synonymous mutation
                            n += 1.0/3.0;
                        }
                    }
                }
            }
            self.ns_sites.insert(codon.clone(), n as f32);
        }
    }

    fn find_mutations(&self,
                      gene: &bio::io::gff::Record,
                      variant_abundances: &Vec<HashMap<String, f32>>,
                      ref_sequence: &Vec<u8>,
                      depth: &Vec<usize>) {
        let strand = gene.strand().expect("No strandedness found");

        // bio::gff documentation says start and end positions are 1-based, so we minus 1
        // Additionally, end position is non-inclusive
        let start = gene.start().clone() as usize - 1;
        let end = gene.end().clone() as usize - 1;
        let frame: usize = gene.frame().parse().unwrap();
        let gene_sequence = ref_sequence[start..end].to_vec();
        debug!("Gene Seq {:?}", String::from_utf8_lossy(&gene_sequence));
        let codon_sequence = get_codons(gene_sequence, frame, strand);

        // Calculate N and S
        let mut N: f32 = 0.0;
        let mut S: f32 = 0.0;
        for codon in codon_sequence {
            let n = self.ns_sites[&codon];
            N += n;
            S += 3.0 - n;
        }
        // Create Nd and Sd values
        let mut Nd: f32 = 0.0;
        let mut Sd: f32 = 0.0;
        // To calculate dN/dS values we need to follow Morelli et al. 2013 and their adaption
        // for dN/dS calculations when using NGS reads also outlined here:
        // http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
        let mut codon = vec!();
        let mut new_codon = vec!();
        for (gene_cursor, variant_map) in variant_abundances[start..end].to_vec().iter().enumerate() {
            let codon_idx = gene_cursor / 3 as usize;
            let codon_cursor = gene_cursor % 3;

            if codon_cursor == 0 {
                if (codon.len() == 3) & (new_codon.len()) == 3 {
                    // get indices of different locations
                    let mut pos = 0 as usize;
                    let mut diffs = vec!();
                    for (c1, c2) in codon.iter().zip(new_codon.iter()) {
                        if c1 != c2 {
                            diffs.push(pos);
                        }
                        pos += 1;
                    }
                    // get permuations of positions
                    let heap = Heap::new(&mut diffs);
                    let mut permutations = Vec::new();
                    for data in heap {
                        permutations.push(data.clone());
                    }
                    // calculate synonymous and non-synonymous for each permutation
                    let mut ns = 0;
                    let mut ss = 0;
                    for permutation in permutations {
                        let mut shifting = codon.clone();
                        let mut old_shift;
                        for pos in permutation {
                            old_shift = shifting.clone();
                            shifting[pos] = new_codon[pos];
                            if self.aminos[&old_shift] != self.aminos[&shifting] {
                                ns += 1;
                            } else {
                                ss += 1;
                            }
                        }
                    }
                    let nd = ns as f32 / permutations.len() as f32;
                    let sd = ss as f32 / permutations.len() as f32;
                    Nd += nd;
                    Sd += sd;

                }
                // begin working on new codon
                codon = codon_sequence[codon_idx].clone();
                new_codon = codon.clone();
            }

            if variant_map.len() > 0 {
                for variant in variant_map.keys() {
                    if variant.len() > 1 {
                        debug!("Frameshift mutation variant {:?}", variant);
                        continue
                    }

                    let variant = variant.as_bytes().to_vec();
                    let codon = codon_sequence[codon_idx].clone();
                    let mut new_codon = codon.clone();
                    new_codon[codon_cursor] = variant[0];
                }
            }
        }
    }
}

pub fn get_codons(sequence: Vec<u8>, frame: usize, strandedness: strand::Strand) -> Vec<Vec<u8>> {

    let codons = match strandedness{
        strand::Strand::Forward | strand::Strand::Unknown => {
            sequence[0+frame..].chunks(3).map(|chunk| chunk.to_vec()).collect::<Vec<Vec<u8>>>()
        },
        strand::Strand::Reverse => {
            let rc = dna::revcomp(sequence);
            rc[0+frame..].chunks(3).map(|chunk| chunk.to_vec()).collect::<Vec<Vec<u8>>>()
        }
    };
    return codons
}