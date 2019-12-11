use std::collections::{HashMap, HashSet, BTreeMap, BTreeSet};
use itertools::izip;
use bio::alphabets::dna;
use bio_types::strand;
use permutohedron::{Heap};
use bio::io::gff;


pub struct CodonTable {
    aminos: HashMap<Vec<u8>, char>,
    starts: HashMap<Vec<u8>, char>,
    ns_sites: HashMap<Vec<u8>, f32>,
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
            1 => {
                NCBITable {
                    aas: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".to_owned(),
                    starts: "---M------**--*----M---------------M----------------------------".to_owned(),
                    base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".to_owned(),
                    base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".to_owned(),
                    base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".to_owned(),
                }
            },
            11 => {
                NCBITable {
                    aas: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".to_owned(),
                    starts: "---M------**--*----M------------MMMM---------------M------------".to_owned(),
                    base1: "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".to_owned(),
                    base2: "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".to_owned(),
                    base3: "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".to_owned(),
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
                      variant_abundances: &HashMap<i32, BTreeMap<String, (f64, f64)>>,
                      indels: &HashMap<i32, BTreeMap<String, BTreeSet<i64>>>,
                      ref_sequence: &Vec<u8>,
                      depth: &Vec<f64>) -> f32;
}

impl Translations for CodonTable {
    fn get_codon_table(&mut self, table_id: usize) {
        // Convert NCBI format to something more usable
        let ncbi_format = NCBITable::get_translation_table(table_id);
        let mut amino_hash = HashMap::new();
        let mut start_hash = HashMap::new();
        for (aa, s, b1, b2, b3) in izip!(ncbi_format.aas.as_bytes(), ncbi_format.starts.as_bytes(),
                                        ncbi_format.base1.as_bytes(), ncbi_format.base2.as_bytes(),
                                        ncbi_format.base3.as_bytes()) {
            let codon = vec!(*b1, *b2, *b3);
            amino_hash.insert(codon.clone(), *aa as char);
            start_hash.insert(codon, *s as char);
        }
        debug!("Amino hash {:?}", amino_hash);
        self.aminos = amino_hash;
        self.starts = start_hash;

        let nucleotides: Vec<u8> = vec!(65, 84, 67, 71);
        debug!("nucleotides: {} ", String::from_utf8_lossy(&nucleotides));
        for (codon, _aa) in self.aminos.iter(){
            let mut n = 0.0;
            for (pos, cod) in codon.iter().enumerate() {
                for nuc in nucleotides.iter() {
                    if cod == nuc {
                        // Ignore this nucleotide
                        continue
                    } else {
                        let mut codon_shift = codon.clone();
                        // Change one position of the codon
                        codon_shift[pos] = nuc.clone();

                        if self.aminos[codon] != self.aminos[&codon_shift] {
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
                      variant_abundances: &HashMap<i32, BTreeMap<String, (f64, f64)>>,
                      indels: &HashMap<i32, BTreeMap<String, BTreeSet<i64>>>,
                      ref_sequence: &Vec<u8>,
                      _depth: &Vec<f64>) -> f32 {
        let strand = gene.strand().expect("No strandedness found");

        // bio::gff documentation says start and end positions are 1-based, so we minus 1
        // Additionally, end position is non-inclusive
        let start = gene.start().clone() as usize - 1;
        let end = gene.end().clone() as usize;
        let frame: usize = gene.frame().parse().unwrap();
        let gene_sequence = ref_sequence[start..end].to_vec();
        debug!("Gene Seq {:?}", String::from_utf8_lossy(&gene_sequence));
        let codon_sequence = get_codons(&gene_sequence, frame, strand);
        debug!("Codon Sequence {:?}", codon_sequence);

        // Calculate N and S
        let mut big_n: f32 = 0.0;
        let mut big_s: f32 = 0.0;
        for codon in codon_sequence.iter() {
            if String::from_utf8(codon.clone()).expect("Unable to interpret codon").contains("N") {
                continue
            } else {
                let n = self.ns_sites[codon];
                big_n += n;
                big_s += 3.0 - n;
            }
        }

        debug!("getting ns_sites N {} S {}", big_n, big_s);

        // Create Nd and Sd values
        let mut big_nd: f32 = 0.0;
        let mut big_sd: f32 = 0.0;

        // dN/dS calculations when using NGS reads outlined here:
        // http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
        // Note, we don't normalize for depth here and instead just use Jukes-Cantor model
        let mut codon: Vec<u8> = vec!();
        let mut new_codons: Vec<Vec<u8>> = vec!();
        let mut positionals = 0;
        let mut total_variants = 0;
        let mut indel_map = BTreeMap::new();
        let dummy = BTreeMap::new();
        for (gene_cursor, cursor) in (start..end).into_iter().enumerate() {
            let variant_map = match variant_abundances.get(&(cursor as i32)){
                Some(map) => map,
                None => &dummy,
            };

            if indels.contains_key(&(cursor as i32)){
                indel_map = indels.get(&(cursor as i32))
                    .expect("No INDEL at this location").to_owned();
            };

            let codon_idx = gene_cursor / 3 as usize;
            let codon_cursor = gene_cursor % 3;
            if String::from_utf8(codon.clone())
                .expect("Unable to interpret codon").contains("N") {
                continue
            }

            if codon_cursor == 0
                || new_codons.len() == 0 {
                for new_codon in new_codons.iter_mut() {
                    if (codon.len() == 3)
                        && (new_codon.len() == 3)
                        && (codon != *new_codon){
                        // get indices of different locations
                        let mut pos = 0 as usize;
                        let mut diffs = vec!();
                        for (c1, c2) in codon.iter().zip(new_codon.iter()) {
                            if c1 != c2 {
                                diffs.push(pos);
                            }
                            pos += 1;
                        }
                        total_variants += diffs.len();
                        // get permuations of positions
                        let heap = Heap::new(&mut diffs);
                        let mut permutations = Vec::new();
                        for data in heap {
                            permutations.push(data.clone());
                        }
                        // calculate synonymous and non-synonymous for each permutation
                        let mut ns = 0;
                        let mut ss = 0;
                        debug!("positional difference {:?} permutations {:?}", diffs, permutations.len());
                        positionals += permutations.len();
                        for permutation in permutations.iter() {
                            let mut shifting = codon.clone();
                            let mut old_shift;
                            for pos in permutation {
                                // Check if one amino acid change causes an syn or non-syn
                                old_shift = shifting.clone();
                                shifting[*pos] = new_codon[*pos];
                                debug!("Old shift {:?}, new {:?}", old_shift, shifting);
                                if self.aminos[&old_shift] != self.aminos[&shifting] {
                                    ns += 1;
                                } else {
                                    ss += 1;
                                }
                            }
                        }
                        let nd = ns as f32 / permutations.len() as f32;
                        let sd = ss as f32 / permutations.len() as f32;
                        big_nd += nd;
                        big_sd += sd;
                    }
                }
                // begin working on new codon
                debug!("Codon idx {} codonds {} gene length {} new_codons {:?}",
                         codon_idx, codon_sequence.len(), gene_sequence.len(), new_codons);
                if codon_sequence.len() == 268 {
                    debug!("{:?}", codon_sequence);
                }
                codon = codon_sequence[codon_idx].clone();
                if codon.len() != 3 {
                    continue
                }
                new_codons = Vec::new();
                new_codons.push(codon.clone());
            }
            if variant_map.len() > 0 {
                debug!("variant map {:?}", variant_map);
                let mut variant_count = 0;
                let mut variant_vec: Vec<_> = variant_map.iter().collect();
                // We look at the most abundant variant first for consistency
                variant_vec.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
                for (variant, _frac) in variant_vec {
                    if indel_map.contains_key(variant)
                        || variant.to_owned().contains("R") {
                        // Frameshift mutations are not included in dN/dS calculations?
                        // Seems weird, but all formulas say no
                        debug!("Frameshift mutation variant {:?}", variant);
                        continue
                    }
                    let variant = variant.as_bytes().to_vec();
                    if variant_count > 0 {
                        // Create a copy of codon up to this point
                        // Not sure if reusing previous variants is bad, but
                        // not doing so can cause randomness to dN/dS values

                        new_codons.push(codon.clone());

                        new_codons[variant_count][codon_cursor] = variant[0];

                        debug!("multi variant codon {:?}", new_codons);
                    } else {

                        for var_idx in 0..new_codons.len() {

                            new_codons[var_idx][codon_cursor] = variant[0];
                        }
                    }
                    variant_count += 1;
                }
            }
        }

        debug!("Nd {} N {}, Sd {} S {} total permutations {} variants {}",
               big_nd, big_n, big_sd, big_s, positionals, total_variants);
        let mut pn = big_nd/big_n;
        let mut ps = big_sd/big_s;
        debug!("pn {} ps {}", pn, ps);
        // Weirdly in the Jukes-Cantor model if pn or ps are 0.75 then the nat log does not resolve
        // No one talks about this in the literature for some reason
        if pn == 0.75 {
            pn = 0.7499
        }
        if ps == 0.75 {
            ps = 0.7499
        }
        let d_n = -(3.0/4.0)*(1.0-(4.0*pn)/3.0).ln();
        let d_s = -(3.0/4.0)*(1.0-(4.0*ps)/3.0).ln();
        debug!("dN {} dS {}", d_n, d_s);
        let mut dnds = d_n/d_s;

        // negative dnds values make no sense, but occur nonetheless
        // Just make them 0.0
        if dnds < 0.0 {
            dnds = 0.0
        }

        return dnds
    }
}

pub fn get_codons(sequence: &Vec<u8>, frame: usize, strandedness: strand::Strand) -> Vec<Vec<u8>> {

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


#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::gff::GffType;

    #[test]
    fn test_floor_division() {
        assert_eq!(0/3 as usize, 0);
        assert_eq!(803/3 as usize, 267);
    }

    #[test]
    fn test_dnds() {
        let mut codon_table = CodonTable::setup();
        codon_table.get_codon_table(11);

        let mut gene_records
            = gff::Reader::from_file("tests/data/dnds.gff", bio::io::gff::GffType::GFF3).expect("Incorrect file path");
        let mut variant_abundances: HashMap<i32, BTreeMap<String, (f64, f64)>> = HashMap::new();
        variant_abundances.insert(13, BTreeMap::new());
        variant_abundances.insert(14, BTreeMap::new());
        let hash = variant_abundances.entry(7).or_insert(BTreeMap::new());
        hash.insert("G".to_string(), (5., 10.));
        hash.insert("R".to_string(), (5., 10.));

        let hash = variant_abundances.entry(11).or_insert(BTreeMap::new());
        hash.insert("C".to_string(), (5., 10.));
        hash.insert("R".to_string(), (5., 10.));

        let hash = variant_abundances.entry(13).or_insert(BTreeMap::new());
        hash.insert("A".to_string(), (5., 10.));
        hash.insert("R".to_string(), (5., 10.));

        let hash = variant_abundances.entry(14).or_insert(BTreeMap::new());
        hash.insert("C".to_string(), (5., 10.));
        hash.insert("R".to_string(), (5., 10.));

        for gene_record in gene_records.records() {
            let gene_record = gene_record.unwrap();
            let ref_sequence = "ATGAAACCCGGGTTTTAA".as_bytes().to_vec();

            let dnds = codon_table.find_mutations(
                &gene_record,
                &variant_abundances,
                &HashMap::new(),
                &ref_sequence,
                &Vec::new());
            assert_eq!(format!("{:.4}", dnds), format!("{}", 0.1247));
        }

    }
}