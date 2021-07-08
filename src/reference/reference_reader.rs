use bio::io::fasta::IndexedReader;
use coverm::genomes_and_contigs::GenomesAndContigs;
use std::fs::File;
use std::collections::HashMap;
use rayon::prelude::*;
use reference::reference_reader_utils::ReferenceReaderUtils;

/**
* Struct handling methods to read and handle information for references
*/
#[derive(Debug)]
pub struct ReferenceReader<'a> {
    indexed_reader: IndexedReader<File>,
    pub current_sequence: Vec<u8>,
    pub genomes_and_contigs: &'a GenomesAndContigs,
    target_names: HashMap<usize, Vec<u8>>,
}

impl<'a> ReferenceReader<'a> {
    pub fn new(
        concatenated_genomes: &Option<String>,
        genomes_and_contigs: &'a GenomesAndContigs,
        number_of_contigs: usize,
    ) -> ReferenceReader<'a> {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(concatenated_genomes);

        ReferenceReader {
            indexed_reader,
            current_sequence: Vec::new(),
            target_names: HashMap::new(),
            genomes_and_contigs: genomes_and_contigs,
        }
    }

    pub fn new_with_target_names(
        concatenated_genomes: &Option<String>,
        genomes_and_contigs: &'a GenomesAndContigs,
        target_names: Vec<&[u8]>,
    ) -> ReferenceReader<'a> {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(concatenated_genomes);

        ReferenceReader {
            indexed_reader,
            current_sequence: Vec::new(),
            target_names: target_names.into_par_iter().enumerate()
                .map(|(tid, target)| (tid, target.to_vec())).collect::<HashMap<usize, Vec<u8>>>(),
            genomes_and_contigs: genomes_and_contigs,
        }
    }

    pub fn add_target(&mut self, target: &[u8], tid: usize) {
        self.target_names.insert(tid, target.to_vec());
    }

    pub fn retrieve_reference_stem(&self, ref_idx: usize) -> String {
        self.genomes_and_contigs.genomes[ref_idx].clone()
    }

    pub fn update_current_sequence_capacity(&mut self, size: usize) {
        self.current_sequence = Vec::with_capacity(size);
    }

    pub fn update_current_sequence_without_capcity(&mut self) {
        self.current_sequence = Vec::new();
    }

    pub fn match_target_name_and_ref_idx(&self, ref_idx: usize, target_name: &String) -> bool {
        match self.genomes_and_contigs.contig_to_genome.get(target_name) {
            Some(ref_id) => *ref_id == ref_idx,
            None => false,
        }
    }

    pub fn retrieve_contig_name_from_tid(&self, tid: usize) -> Option<&Vec<u8>> {
        return self.target_names.get(&tid)
    }

    pub fn fetch_contig_from_reference_by_contig_name(
        &mut self,
        contig_name: &[u8],
        ref_idx: usize,
    ) {
        match self.indexed_reader.fetch_all(std::str::from_utf8(contig_name).unwrap()) {
            Ok(reference) => reference,
            Err(_e) => match self.indexed_reader.fetch_all(&format!(
                "{}~{}",
                self.genomes_and_contigs.genomes[ref_idx],
                std::str::from_utf8(contig_name).unwrap()
            )) {
                Ok(reference) => reference,
                Err(e) => {
                    println!(
                        "Cannot read sequence from reference {} {:?}",
                        format!(
                            "{}~{}",
                            &self.genomes_and_contigs.genomes[ref_idx],
                            std::str::from_utf8(contig_name).unwrap()
                        ),
                        e,
                    );
                    std::process::exit(1);
                }
            },
        };
    }

    pub fn fetch_contig_from_reference_by_tid(
        &mut self,
        tid: usize,
        ref_idx: usize,
    ) {
        match self.indexed_reader.fetch_all(std::str::from_utf8(&self.target_names[&tid]).unwrap()) {
            Ok(reference) => reference,
            Err(_e) => match self.indexed_reader.fetch_all(&format!(
                "{}~{}",
                self.genomes_and_contigs.genomes[ref_idx],
                std::str::from_utf8(&self.target_names[&tid]).unwrap()
            )) {
                Ok(reference) => reference,
                Err(e) => {
                    println!(
                        "Cannot read sequence from reference {} {:?}",
                        format!(
                            "{}~{}",
                            &self.genomes_and_contigs.genomes[ref_idx],
                            std::str::from_utf8(&self.target_names[&tid]).unwrap()
                        ),
                        e,
                    );
                    std::process::exit(1);
                }
            },
        };
    }

    pub fn read_sequence_to_vec(
        &mut self,
    ) {
        match self.indexed_reader.read(&mut self.current_sequence) {
            Ok(reference) => reference,
            Err(e) => {
                println!(
                    "Cannot read sequence from reference {:?}",
                    e,
                );
                std::process::exit(1)
            }
        };
    }

    /**
    * Takes a contig name &[u8] and splits around a suspected separator character
    */
    pub fn split_contig_name<'b>(contig_name: &'b [u8], separator: u8) -> &'b [u8] {
        match contig_name.into_par_iter().position_first(|&x| x == separator) {
            Some(position) => {
                return &contig_name[position..contig_name.len()]
            },
            None => {
                return contig_name
            }
        }
    }
}