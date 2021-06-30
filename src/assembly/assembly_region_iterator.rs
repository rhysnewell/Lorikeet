use utils::reference_reader_utils::ReferenceReader;
use activity_profile::band_pass_activity_profile::BandPassActivityProfile;
use haplotype::haplotype_caller_engine::HaplotypeCallerEngine;
use assembly::assembly_region::AssemblyRegion;
use activity_profile::activity_profile::{ActivityProfile, Profile};
use reads::bird_tool_reads::BirdToolRead;
use coverm::bam_generator::{generate_indexed_named_bam_readers_from_bam_files, IndexedNamedBamReader};
use rust_htslib::bam::Record;
use coverm::FlagFilter;

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
#[derive(Debug, Clone)]
pub struct AssemblyRegionIterator {
    reference_reader: &mut ReferenceReader,
    indexed_bam_readers: &Vec<String>,
    pub pending_regions: Vec<AssemblyRegion>,
    n_threads: u32,
    // previous_regions_reads: Vec<BirdToolRead>,
}

impl AssemblyRegionIterator {
    pub fn new(
        reference_reader: &mut ReferenceReader,
        activity_profile: &mut ActivityProfile,
        indexed_bam_readers: &Vec<String>,
        assembly_region_padding: usize,
        min_assembly_region_size: usize,
        max_assembly_region_size: usize,
        n_threads: u32
    ) -> AssemblyRegionIterator {

        // Assume no forced conversion here since we have already traverse the entire
        // activity profile prior to reaching here. This is quite different to how
        // GATK handles it but I assume it ends up working the same?

        let mut pending_regions = activity_profile.pop_ready_assembly_regions(
            assembly_region_padding,
            min_assembly_region_size,
            max_assembly_region_size,
            false
        );
        AssemblyRegionIterator {
            reference_reader,
            indexed_bam_readers,
            pending_regions,
            n_threads,
        }
    }

    pub fn fill_next_assembly_region_with_reads(
        region: &mut AssemblyRegion,
        flag_filters: &FlagFilter,
        n_threads: u32
    ) {
        // We don't need to check previous region reads as the implementation of fetch we have
        // should retrieve all reads regardles of if they have been seen before
        let mut record = Record::new(); // Empty bam record
        indexed_bam_readers.iter().enumerate().for_each(
            |(sample_idx, bam_generator)| {
                let bam_generated = generate_indexed_named_bam_readers_from_bam_files(
                    vec![&bam_generator],
                    n_threads,
                )
                    .into_iter()
                    .next()
                    .unwrap();

                bam_generated.fetch((
                    region.get_contig() as i32,
                    region.get_start() as i64,
                    region.get_end() as i64,
                ));

                while bam_generated.read(&mut record) == true {
                    // TODO: Implement read filtering here
                    if (!flag_filters.include_supplementary
                        && record.is_supplementary()) // We want supp alignments for longreads
                        || (!flag_filters.include_secondary
                        && record.is_secondary())
                        || record.is_unmapped()
                    // Check against filter flags and current sample type
                    {
                        skipped_reads += 1;
                        continue;
                    } else {
                        region.add(BirdToolRead::new(*record.clone()))
                    }
                }
            });
    }
}

impl Iterator for AssemblyRegionIterator {
    type Item = AssemblyRegion;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.pending_regions.next())
    }
}