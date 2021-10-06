use assembly::assembly_region::AssemblyRegion;
use coverm::bam_generator::{
    generate_indexed_named_bam_readers_from_bam_files, IndexedNamedBamReader,
};
use coverm::FlagFilter;
use estimation::lorikeet_engine::ReadType;
use rayon::prelude::*;
use reads::bird_tool_reads::BirdToolRead;
use rust_htslib::bam::Record;

/**
 * Given a {@link BandPassActivityProfile} and {@link AssemblyRegionWalker}, iterates over each {@link AssemblyRegion} within
 * that shard, using the provided {@link AssemblyRegionEvaluator} to determine the boundaries between assembly
 * regions.
 *
 * Loads the reads from the shard as lazily as possible to minimize memory usage.
 *
 * This iterator represents the core of the {@link AssemblyRegionWalker} traversal.
 *
 * NOTE: the provided shard must have appropriate read filters set on it for this traversal type (ie., unmapped
 * and malformed reads must be filtered out).
 */
#[derive(Debug)]
pub struct AssemblyRegionIterator<'a> {
    indexed_bam_readers: &'a Vec<String>,
    n_threads: u32,
    // previous_regions_reads: Vec<BirdToolRead>,
}

impl<'a> AssemblyRegionIterator<'a> {
    pub fn new(indexed_bam_readers: &'a Vec<String>, n_threads: u32) -> AssemblyRegionIterator<'a> {
        // Assume no forced conversion here since we have already traverse the entire
        // activity profile prior to reaching here. This is quite different to how
        // GATK handles it but I assume it ends up working the same?
        AssemblyRegionIterator {
            indexed_bam_readers,
            n_threads,
        }
    }

    pub fn fill_next_assembly_region_with_reads(
        &self,
        region: &mut AssemblyRegion,
        flag_filters: &FlagFilter,
        n_threads: u32,
        short_read_bam_count: usize,
        long_read_bam_count: usize,
    ) {
        // We don't need to check previous region reads as the implementation of fetch we have
        // should retrieve all reads regardles of if they have been seen before
        let mut records: Vec<BirdToolRead> = self
            .indexed_bam_readers
            .par_iter()
            .enumerate()
            .flat_map(|(sample_idx, bam_generator)| {
                let mut record = Record::new(); // Empty bam record
                let mut bam_generated = generate_indexed_named_bam_readers_from_bam_files(
                    vec![&bam_generator],
                    n_threads,
                )
                .into_iter()
                .next()
                .unwrap();

                bam_generated.fetch((
                    region.get_contig() as i32,
                    region.get_padded_span().start as i64, // fetch is possibly 1-based
                    region.get_padded_span().end as i64,
                ));

                let read_type = if sample_idx < short_read_bam_count {
                    ReadType::Short
                } else {
                    ReadType::Long
                };

                let mut records = Vec::new(); // container for the records to be collected

                while bam_generated.read(&mut record) == true {
                    // TODO: Implement read filtering here
                    if (!flag_filters.include_supplementary
                        && record.is_supplementary()
                        && read_type != ReadType::Long) // We want supp alignments for longreads
                        || (!flag_filters.include_secondary
                        && record.is_secondary())
                        || record.is_unmapped()
                    // Check against filter flags and current sample type
                    {
                        continue;
                    } else {
                        records.push(BirdToolRead::new(record.clone(), sample_idx, read_type));
                    };
                }

                records
            })
            .collect::<Vec<BirdToolRead>>();
        debug!("Adding {} reads to region", records.len());
        records.par_sort_unstable();
        region.add_all(records);
    }
}

// impl Iterator for AssemblyRegionIterator<'_> {
//     type Item = AssemblyRegion;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         Some(self.pending_regions.next())
//     }
// }
