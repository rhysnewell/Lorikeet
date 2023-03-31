use bio::io::fasta::IndexedReader;
use hashlink::LinkedHashSet;
use std::cmp::min;
use std::collections::HashMap;
use std::fs::File;

use crate::reference::reference_reader_utils::GenomesAndContigs;
use crate::reference::reference_reader_utils::ReferenceReaderUtils;
use crate::utils::simple_interval::{Locatable, SimpleInterval};

/**
* Struct handling methods to read and handle information for references
* indexed_reader represents the #[IndexedReader<File>]
* current_sequence is the byte array current being read from the indexed_reader file
* genomes_and_contigs holds the #[GenomesAndContigs] struct matching contig string names to tids
* reference index to tid holds a reference index matched to a LinkedHashSet of all associated tids
* target_names matches the tid to a byte array representation of the contigs name in the Fasta file
* target_lens matches the tid of a contig to its observed length in the Fasta file in base pairs
*/
#[derive(Debug)]
pub struct ReferenceReader {
    indexed_reader: IndexedReader<File>,
    pub current_sequence: Vec<u8>,
    pub genomes_and_contigs: GenomesAndContigs,
    reference_index_to_tid: HashMap<usize, LinkedHashSet<usize>>,
    target_names: HashMap<usize, Vec<u8>>,
    pub target_lens: HashMap<usize, u64>,
    genome_path: Option<String>,
}

impl Clone for ReferenceReader {
    fn clone(&self) -> Self {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(&self.genome_path);

        Self {
            indexed_reader,
            // current_sequence: self.current_sequence.clone(),
            current_sequence: Vec::new(),
            genomes_and_contigs: self.genomes_and_contigs.clone(),
            reference_index_to_tid: self.reference_index_to_tid.clone(),
            target_names: self.target_names.clone(),
            target_lens: self.target_lens.clone(),
            genome_path: self.genome_path.clone(),
        }
    }
}

impl ReferenceReader {
    pub fn new(
        concatenated_genomes: &Option<String>,
        genomes_and_contigs: GenomesAndContigs,
        _number_of_contigs: usize,
    ) -> ReferenceReader {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(concatenated_genomes);

        ReferenceReader {
            indexed_reader,
            current_sequence: Vec::new(),
            reference_index_to_tid: HashMap::new(),
            target_names: HashMap::new(),
            target_lens: HashMap::new(),
            genomes_and_contigs,
            genome_path: concatenated_genomes.clone(),
        }
    }

    pub fn new_with_target_names(
        concatenated_genomes: &Option<String>,
        genomes_and_contigs: GenomesAndContigs,
        target_names: Vec<&[u8]>,
    ) -> Self {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(concatenated_genomes);

        Self {
            indexed_reader,
            current_sequence: Vec::new(),
            reference_index_to_tid: HashMap::new(),
            target_names: target_names
                .into_iter()
                .enumerate()
                .map(|(tid, target)| (tid, target.to_vec()))
                .collect::<HashMap<usize, Vec<u8>>>(),
            genomes_and_contigs,
            target_lens: HashMap::new(),
            genome_path: concatenated_genomes.clone(),
        }
    }

    /// Somewhat efficient cloning function. Only takes the info needed for a given ref idx
    pub fn new_from_reader_with_reference_index(reader: &Self, ref_idx: usize) -> Self {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(&reader.genome_path);

        if reader.reference_index_to_tid.contains_key(&ref_idx) {
            let ref_tids = reader.reference_index_to_tid.get(&ref_idx).unwrap().clone();
            let mut target_names = HashMap::with_capacity(ref_tids.len());
            let mut target_lens = HashMap::with_capacity(ref_tids.len());
            ref_tids.iter().for_each(|tid| {
                target_names.insert(*tid, reader.target_names.get(tid).unwrap().clone());
                target_lens.insert(*tid, reader.target_lens.get(tid).unwrap().clone());
            });

            let mut ref_idx_to_tid = HashMap::new();
            ref_idx_to_tid.insert(ref_idx, ref_tids);

            let genomes_and_contigs = GenomesAndContigs::new();

            Self {
                indexed_reader,
                current_sequence: Vec::new(),
                reference_index_to_tid: ref_idx_to_tid,
                genomes_and_contigs,
                target_lens,
                target_names,
                genome_path: reader.genome_path.clone(),
            }
        } else {
            reader.clone()
        }
    }

    pub fn new_from_reader_with_tid_and_rid(
        reader: &ReferenceReader,
        ref_idx: usize,
        tid: usize,
    ) -> Self {
        let indexed_reader = ReferenceReaderUtils::retrieve_reference(&reader.genome_path);

        if reader.reference_index_to_tid.contains_key(&ref_idx) {
            let mut target_names = HashMap::with_capacity(1);
            let mut target_lens = HashMap::with_capacity(1);
            target_names.insert(tid, reader.target_names.get(&tid).unwrap().clone());
            target_lens.insert(tid, reader.target_lens.get(&tid).unwrap().clone());

            let mut ref_idx_to_tid = HashMap::new();
            let mut tids = LinkedHashSet::with_capacity(1);
            tids.insert(tid);
            ref_idx_to_tid.insert(ref_idx, tids);

            let genomes_and_contigs = GenomesAndContigs::new();

            Self {
                indexed_reader,
                current_sequence: Vec::new(),
                reference_index_to_tid: ref_idx_to_tid,
                genomes_and_contigs,
                target_lens,
                target_names,
                genome_path: reader.genome_path.clone(),
            }
        } else {
            reader.clone()
        }
    }

    pub fn get_target_name(&self, tid: usize) -> &[u8] {
        self.target_names.get(&tid).unwrap()
    }

    pub fn add_target(&mut self, target: &[u8], tid: usize) {
        self.target_names.insert(tid, target.to_vec());
    }

    pub fn update_ref_index_tids(&mut self, ref_index: usize, tid: usize) {
        let tids = self
            .reference_index_to_tid
            .entry(ref_index)
            .or_insert_with(LinkedHashSet::new);
        tids.insert(tid);
    }

    pub fn retrieve_tids_for_ref_index(&self, ref_index: usize) -> Option<&LinkedHashSet<usize>> {
        self.reference_index_to_tid.get(&ref_index)
    }

    pub fn add_length(&mut self, tid: usize, length: u64) {
        self.target_lens.insert(tid, length);
    }

    pub fn add_lengths(&mut self, target_lengths: HashMap<usize, u64>) {
        self.target_lens = target_lengths
    }

    pub fn get_contig_length(&self, tid: usize) -> u64 {
        *self.target_lens.get(&tid).unwrap_or(&0)
    }

    pub fn retrieve_reference_stem(&self, ref_idx: usize) -> String {
        self.genomes_and_contigs.genomes[ref_idx].clone()
    }

    pub fn update_current_sequence_capacity(&mut self, size: usize) {
        self.current_sequence = Vec::with_capacity(size);
    }

    pub fn update_current_sequence_without_capacity(&mut self) {
        self.current_sequence = Vec::new();
    }

    pub fn retrieve_contig_name_from_tid(&self, tid: usize) -> Option<&Vec<u8>> {
        return self.target_names.get(&tid);
    }

    pub fn fetch_contig_from_reference_by_contig_name(
        &mut self,
        contig_name: &[u8],
        ref_idx: usize,
    ) {
        match self
            .indexed_reader
            .fetch_all(std::str::from_utf8(contig_name).unwrap())
        {
            Ok(reference) => reference,
            Err(_e) => match self.indexed_reader.fetch_all(&format!(
                "{}",
                std::str::from_utf8(contig_name)
                    .unwrap()
                    .splitn(2, '~')
                    .nth(1)
                    .unwrap()
            )) {
                Ok(reference) => reference,
                Err(_e) => {
                    match self.indexed_reader.fetch_all(&format!(
                        "{}~{}",
                        self.genomes_and_contigs.genomes[ref_idx],
                        std::str::from_utf8(contig_name).unwrap()
                    )) {
                        Ok(reference) => reference,
                        Err(e) => {
                            panic!(
                                "Cannot read sequence from {}: {}",
                                std::str::from_utf8(contig_name).unwrap(),
                                e
                            )
                        }
                    }
                }
            },
        };
    }

    pub fn fetch_contig_from_reference_by_tid(
        &mut self,
        tid: usize,
        ref_idx: usize,
    ) -> Result<(), std::io::Error> {
        match self
            .indexed_reader
            .fetch_all(std::str::from_utf8(&self.target_names[&tid]).unwrap())
        {
            Ok(reference) => Ok(reference),
            Err(_e) => match self.indexed_reader.fetch_all(&format!(
                "{}",
                match std::str::from_utf8(&self.target_names[&tid])
                    .unwrap()
                    .splitn(2, '~')
                    .nth(1)
                {
                    None => std::str::from_utf8(&self.target_names[&tid]).unwrap(),
                    Some(contig) => contig,
                }
            )) {
                Ok(reference) => Ok(reference),
                Err(_e) => {
                    match self.indexed_reader.fetch_all(&format!(
                        "{}~{}",
                        self.genomes_and_contigs.genomes[ref_idx],
                        std::str::from_utf8(&self.target_names[&tid]).unwrap()
                    )) {
                        Ok(reference) => Ok(reference),
                        Err(e) => {
                            debug!(
                                "Cannot read sequence from {}: {}",
                                std::str::from_utf8(&self.target_names[&tid]).unwrap(),
                                e
                            );
                            Err(e)
                        }
                    }
                }
            },
        }
    }

    /// Fetches the reference sequence from a given SimpleInterval
    /// The return position is 0-base start and stop inclusive
    pub fn fetch_reference_context(&mut self, ref_idx: usize, interval: &SimpleInterval) {
        match self.indexed_reader.fetch(
            std::str::from_utf8(&self.target_names[&interval.get_contig()]).unwrap(),
            interval.get_start() as u64,
            min(
                (interval.get_end() + 1) as u64,
                *self.target_lens.get(&interval.get_contig()).unwrap(),
            ),
        ) {
            Ok(reference) => reference,
            Err(_e) => match self.indexed_reader.fetch(
                &format!(
                    "{}",
                    std::str::from_utf8(&self.target_names[&interval.get_contig()])
                        .unwrap()
                        .splitn(2, '~')
                        .nth(1)
                        .unwrap()
                ),
                interval.get_start() as u64,
                min(
                    (interval.get_end() + 1) as u64,
                    *self.target_lens.get(&interval.get_contig()).unwrap(),
                ),
            ) {
                Ok(reference) => reference,
                Err(_e) => match self.indexed_reader.fetch(
                    &format!(
                        "{}~{}",
                        self.genomes_and_contigs.genomes[ref_idx],
                        std::str::from_utf8(&self.target_names[&interval.get_contig()]).unwrap()
                    ),
                    interval.get_start() as u64,
                    min(
                        (interval.get_end() + 1) as u64,
                        *self.target_lens.get(&interval.get_contig()).unwrap(),
                    ),
                ) {
                    Ok(reference) => reference,
                    Err(e) => {
                        panic!(
                            "Cannot read sequence from reference {} {:?}",
                            format!(
                                "{}~{}",
                                &self.genomes_and_contigs.genomes[ref_idx],
                                std::str::from_utf8(&self.target_names[&interval.get_contig()])
                                    .unwrap()
                            ),
                            e,
                        );
                    }
                },
            },
        };
    }

    pub fn read_sequence_to_vec(&mut self) {
        match self.indexed_reader.read(&mut self.current_sequence) {
            Ok(reference) => reference,
            Err(e) => {
                panic!("Cannot read sequence from reference {:?}", e,);
            }
        };
    }

    /**
     * Takes a contig name &[u8] and splits around a suspected separator character
     */
    pub fn split_contig_name<'b>(contig_name: &'b [u8], separator: u8) -> &'b [u8] {
        match contig_name.into_iter().position(|&x| x == separator) {
            Some(position) => return &contig_name[(position + 1)..contig_name.len()], // + 1 because we do  not want to include the sep
            None => return contig_name,
        }
    }
}
