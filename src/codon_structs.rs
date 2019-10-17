use std::collections::{HashMap, HashSet};
use std::str;
use itertools::izip;

pub struct CodonTable {
    aminos: HashMap<char, HashSet<String>>,
    starts: HashMap<char, HashSet<String>>,
}


pub struct NCBITable {
    AAs: String,
    Starts: String,
    Base1: String,
    Base2: String,
    Base3: String,
}


impl NCBITable {
    fn get_translation_table(table_id: usize) -> NCBITable {
        let ret = match table_id {
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
        };
        return ret
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
}

impl Translations for CodonTable {
    fn get_codon_table(&mut self, table_id: usize) {
        let ncbi_format = NCBITable::get_translation_table(table_id);
        let mut amino_hash = HashMap::new();
        let mut start_hash = HashMap::new();
        for (aa, s, b1, b2, b3) in izip!(ncbi_format.AAs.as_bytes(), ncbi_format.Starts.as_bytes(),
                                        ncbi_format.Base1.as_bytes(), ncbi_format.Base2.as_bytes(),
                                        ncbi_format.Base3.as_bytes()) {
            let mut codon = b1.to_string() + str::from_utf8(&[*b2]).expect("Base cannot be read");
            codon = codon + str::from_utf8(&[*b3]).expect("Base cannot be read");
            let amino = amino_hash.entry(*aa as char).or_insert(HashSet::new());
            amino.insert(codon.clone());
            let start = start_hash.entry(*s as char).or_insert(HashSet::new());
            start.insert(codon);
        }
        self.aminos = amino_hash;
        self.starts = start_hash;
    }
}