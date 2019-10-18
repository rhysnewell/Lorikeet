use std::collections::{HashMap, HashSet};
use std::str;
use itertools::{izip, Itertools};
use bio::alphabets::dna;
use bio_types::strand;


pub struct CodonTable {
    aminos: HashMap<Vec<u8>, char>,
    starts: HashMap<Vec<u8>, char>,
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
        }
    }
}

pub trait Translations {
    fn get_codon_table(&mut self, table_id: usize);
    fn find_mutations(&self,
                      gene: &bio::io::gff::Record,
                      variant_abundances: &Vec<HashMap<String, f32>>,
                      ref_sequence: &Vec<u8>);
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
    }

    fn find_mutations(&self,
                      gene: &bio::io::gff::Record,
                      variant_abundances: &Vec<HashMap<String, f32>>,
                      ref_sequence: &Vec<u8>) {
        let strand = gene.strand().expect("No strandedness found");

        // bio::gff documentation says start and end positions are 1-based, so we minus 1
        // Additionally, end position is non-inclusive
        let start = gene.start().clone() as usize - 1;
        let end = gene.end().clone() as usize - 1;
        let frame: usize = gene.frame().parse().unwrap();
        let gene_sequence = ref_sequence[start..end].to_vec();
        debug!("Gene Seq {:?}", String::from_utf8_lossy(&gene_sequence));
        let codon_sequence = get_codons(gene_sequence, frame, strand);
        for (gene_cursor, variant_map) in variant_abundances[start..end].to_vec().iter().enumerate() {
            let codon_idx = gene_cursor / 3 as usize;
            let codon_cursor = gene_cursor % 3;
            if variant_map.len() > 0 {
                for variant in variant_map.keys() {
                    let variant = variant.as_bytes().to_vec();
                    let codon = codon_sequence[codon_idx].clone();
                    let mut new_codon = codon.clone();
                    new_codon[codon_cursor] = variant[0];
                    if variant.len() > 1 {
                        println!("Frameshift mutation variant {:?}", variant);
                    } else if self.aminos[&codon] != self.aminos[&new_codon] {
                        println!("Non-synonymous mutation old amino {} new amino {}",
                                 self.aminos[&codon], self.aminos[&new_codon]);
                    } else if self.aminos[&codon] == self.aminos[&new_codon] {
                        println!("Synonymous mutation old amino {} new amino {}",
                                 self.aminos[&codon], self.aminos[&new_codon]);
                    }
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