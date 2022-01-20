use model::byte_array_allele::ByteArrayAllele;
use model::variant_context::{VariantContext, VariantType};
use reference::reference_reader::ReferenceReader;
use std::collections::{BTreeMap, BinaryHeap};
use std::fs::{create_dir_all, File};
use std::io::Write;
use std::path::Path;
use utils::simple_interval::Locatable;

/// Struct housing methods for writing out genomes when given specific variant information
/// Basically a wrapper for reference reader
pub struct ReferenceWriter<'a> {
    reference_reader: ReferenceReader,
    output_prefix: &'a str,
}

impl<'a> ReferenceWriter<'a> {
    pub fn new(reference_reader: ReferenceReader, output_prefix: &'a str) -> ReferenceWriter<'a> {
        // ensure path exists
        create_dir_all(output_prefix).expect("Unable to create output directory");
        Self {
            reference_reader,
            output_prefix,
        }
    }

    /// Generates the potential strain genomes calculated by Lorikeet. The VariantContexts are expected
    /// To be tagged with one or more strain genomes in their `attributes` with `VariantAnnotation::Strain`
    /// tag
    pub fn generate_strains(
        &mut self,
        variant_contexts: Vec<VariantContext>,
        ref_idx: usize,
        strain_ids_present: Vec<usize>,
    ) {
        let mut grouped_variant_contexts = Self::split_variant_contexts_by_tid(variant_contexts);
        let mut tids = self
            .reference_reader
            .retrieve_tids_for_ref_index(ref_idx)
            .unwrap()
            .clone();

        for strain_idx in strain_ids_present {
            let file_name = format!(
                "{}/{}_strain_{}.fna",
                self.output_prefix,
                self.reference_reader.genomes_and_contigs.genomes[ref_idx],
                strain_idx,
            );

            let file_path = Path::new(&file_name);
            debug!("File path {}", &file_name);
            // Open new reference file or create one
            let mut file_open =
                File::create(file_path).expect("No Read or Write Permission in current directory");
            for tid in tids.iter() {
                self.reference_reader
                    .fetch_contig_from_reference_by_tid(*tid, ref_idx);
                self.reference_reader.read_sequence_to_vec();

                let mut new_bases = std::mem::take(&mut self.reference_reader.current_sequence);
                let old_length = new_bases.len();
                // This value holds how far right or left the vc location has shifted as we add indels
                let mut offset = 0;
                let mut variant_contexts_of_contig = grouped_variant_contexts.get_mut(&tid);
                let mut variations = 0;
                match variant_contexts_of_contig {
                    Some(variant_contexts_of_contig) => {
                        for mut vc in variant_contexts_of_contig.iter_mut() {
                            if vc.part_of_strain(strain_idx) {
                                let alternate_allele = vc.get_alternate_alleles()[0].clone();
                                let variant_type = vc.get_type().clone();
                                let is_ref = alternate_allele.is_ref;
                                Self::modify_reference_bases_based_on_variant_type(
                                    &mut new_bases,
                                    alternate_allele,
                                    &mut vc,
                                    variant_type,
                                    &mut offset,
                                );
                                variations += if is_ref { 0 } else { 1 };
                            }
                        }
                    }
                    None => {
                        // pass
                    }
                }

                // write the contig header
                writeln!(
                    file_open,
                    ">{} strain_id={} old_length={} new_length={} variations={}",
                    std::str::from_utf8(self.reference_reader.get_target_name(*tid)).unwrap(),
                    strain_idx,
                    old_length,
                    new_bases.len(),
                    variations
                )
                .expect("Unable to write to file");

                // write out the actual contig
                for line in new_bases[..].chunks(60).into_iter() {
                    file_open.write_all(line).unwrap();
                    file_open.write_all(b"\n").unwrap();
                }
            }
        }
    }

    /// Generates the per sample consensus genomes based on the provided variant contexts.
    /// The consensus is defined as the most dominant variant at a given position on the reference
    /// genome measured by read depth.
    pub fn generate_consensus(
        &mut self,
        variant_contexts: Vec<VariantContext>,
        ref_idx: usize,
        samples: &[&str],
    ) {
        let mut grouped_variant_contexts = Self::split_variant_contexts_by_tid(variant_contexts);
        let mut tids = self
            .reference_reader
            .retrieve_tids_for_ref_index(ref_idx)
            .unwrap()
            .clone();
        for (sample_index, sample_name) in samples.into_iter().enumerate() {
            let file_name = format!(
                "{}/{}_consensus_{}.fna",
                self.output_prefix,
                self.reference_reader.genomes_and_contigs.genomes[ref_idx],
                &sample_name.rsplitn(2, '/').next().unwrap(),
            );
            let file_path = Path::new(&file_name);
            debug!("File path {}", &file_name);
            // Open new reference file or create one
            let mut file_open = File::create(file_path).unwrap_or_else(|_| {
                panic!(
                    "No Read or Write Permission in current directory: {:?}",
                    file_path
                )
            });
            for tid in tids.iter() {
                self.reference_reader
                    .fetch_contig_from_reference_by_tid(*tid, ref_idx);
                self.reference_reader.read_sequence_to_vec();
                debug!(
                    "Fetched length {} tid {} ref_idx {} ",
                    self.reference_reader.current_sequence.len(),
                    *tid,
                    ref_idx
                );
                let mut new_bases = std::mem::take(&mut self.reference_reader.current_sequence);
                let old_length = new_bases.len();
                debug!("Contig length {}", old_length);
                // This value holds how far right or left the vc location has shifted as we add indels
                let mut offset = 0;
                let mut variant_contexts_of_contig = grouped_variant_contexts.get_mut(&tid);
                let mut variations = 0;
                match variant_contexts_of_contig {
                    Some(variant_contexts_of_contig) => {
                        for mut vc in variant_contexts_of_contig.iter_mut() {
                            let consensus_allele = vc.get_consensus_allele(sample_index);
                            match consensus_allele {
                                Some(consensus_allele) => {
                                    let variant_type = vc.get_type().clone();
                                    let is_ref = consensus_allele.is_ref;
                                    Self::modify_reference_bases_based_on_variant_type(
                                        &mut new_bases,
                                        consensus_allele,
                                        &mut vc,
                                        variant_type,
                                        &mut offset,
                                    );
                                    variations += if is_ref { 0 } else { 1 };
                                }
                                None => continue,
                            }
                        }
                    }
                    None => {
                        // pass
                    }
                }

                debug!(
                    "Writing contig {}",
                    std::str::from_utf8(self.reference_reader.get_target_name(*tid)).unwrap()
                );
                // write the contig header
                writeln!(
                    file_open,
                    ">{} sample_consensus={} old_length={} new_length={} variations={}",
                    std::str::from_utf8(self.reference_reader.get_target_name(*tid)).unwrap(),
                    sample_name,
                    old_length,
                    new_bases.len(),
                    variations
                )
                .expect("Unable to write to file");

                // write out the actual contig
                for line in new_bases[..].chunks(60).into_iter() {
                    file_open.write_all(line).unwrap();
                    file_open.write_all(b"\n").unwrap();
                }
            }
        }
    }

    /// Takes a list of variant contexts and returns a BTreeMap with contexts grouped by which
    /// contig they appear on. Additionally, the Vector of contexts for each contig is coordinate
    /// sorted so the contexts appear in order in which they occur on the contig
    fn split_variant_contexts_by_tid(
        variant_contexts: Vec<VariantContext>,
    ) -> BTreeMap<usize, Vec<VariantContext>> {
        let mut grouped_variant_contexts = BTreeMap::new();
        for vc in variant_contexts {
            let vcs_on_contig = grouped_variant_contexts
                .entry(vc.loc.get_contig())
                .or_insert_with(BinaryHeap::new);
            vcs_on_contig.push(vc);
        }

        grouped_variant_contexts
            .into_iter()
            .map(|(tid, heap)| (tid, heap.into_sorted_vec()))
            .collect()
    }

    pub fn modify_reference_bases_based_on_variant_type(
        new_bases: &mut Vec<u8>,
        consensus_allele: ByteArrayAllele,
        vc: &mut VariantContext,
        variant_type: VariantType,
        offset: &mut i64,
    ) {
        match variant_type {
            VariantType::Symbolic => {
                if consensus_allele.is_span_del() {
                    // delete from reference, replace with nothing
                    new_bases.splice(
                        ((vc.loc.start as i64 + 1 + *offset) as usize)
                            ..=((vc.loc.end as i64 + *offset) as usize),
                        (0..0).into_iter(),
                    );
                    *offset -= vc.loc.get_length_on_reference() as i64 - 1;
                };
            }
            VariantType::Snp => {
                new_bases[((vc.loc.start as i64 + *offset) as usize)] = consensus_allele.bases[0];
            }
            VariantType::Indel => {
                let allele_len = consensus_allele.bases.len();
                new_bases.splice(
                    ((vc.loc.start as i64 + *offset) as usize)
                        ..=(((vc.loc.start + vc.get_reference().bases.len() - 1) as i64 + *offset)
                            as usize),
                    consensus_allele.bases.into_iter(),
                );

                if vc.loc.get_length_on_reference() == 1 {
                    // insertion
                    *offset += allele_len as i64 - 1;
                } else {
                    *offset -= vc.loc.get_length_on_reference() as i64 - 1;
                };
            }
            VariantType::Mnp => {
                let allele_len = consensus_allele.bases.len();
                new_bases.splice(
                    ((vc.loc.start as i64 + *offset) as usize)
                        ..=(((vc.loc.start + vc.get_reference().bases.len() - 1) as i64 + *offset)
                            as usize),
                    consensus_allele.bases.into_iter(),
                );

                if vc.loc.get_length_on_reference() < allele_len {
                    // gaining bases so increase offset
                    *offset += allele_len as i64 - 1 - vc.loc.get_length_on_reference() as i64;
                } else {
                    // losing bases so decrease offset
                    *offset -= vc.loc.get_length_on_reference() as i64 - allele_len as i64 - 1;
                }
            }
            VariantType::Mixed => {
                // need to determine the type the actual allele came out as
                let new_variant_type = VariantContext::type_of_biallelic_variant(
                    vc.get_reference(),
                    &consensus_allele,
                );
                return Self::modify_reference_bases_based_on_variant_type(
                    new_bases,
                    consensus_allele,
                    vc,
                    new_variant_type,
                    offset,
                );
            }
            _ => {
                // pass on everything else for now.
            }
        }
    }
}
