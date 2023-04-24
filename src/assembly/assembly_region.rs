use rayon::prelude::*;
use std::cmp::min;

use crate::reads::bird_tool_reads::BirdToolRead;
use crate::reads::read_clipper::ReadClipper;
use crate::reference::reference_reader::ReferenceReader;
use crate::utils::interval_utils::IntervalUtils;
use crate::utils::simple_interval::Locatable;
use crate::utils::simple_interval::SimpleInterval;

const MINIMUM_ACTIVITY_DENSITY_THRESHOLD: f32 = 0.2;
const DEFAULT_ADDITIONAL_KMERS: [usize; 3] = [19, 35, 47];

/**
 * Region of the genome that gets assembled by the local assembly engine.
 *
 * As AssemblyRegion is defined by two intervals -- a primary interval containing a territory for variant calling and a second,
 * padded, interval for assembly -- as well as the reads overlapping the padded interval.  Although we do not call variants in the padded interval,
 * assembling over a larger territory improves calls in the primary territory.
 *
 * This concept is complicated somewhat by the fact that these intervals are mutable and the fact that the AssemblyRegion onject lives on after
 * assembly during local realignment during PairHMM.  Here is an example of the life cycle of an AssemblyRegion:
 *
 * Suppose that the HaplotypeCaller engine finds an evidence for a het in a pileup at locus 400 -- that is, it produces
 * an {@code ActivityProfileState} with non-zero probability at site 400 and passes it to its {@code ActivityProfile}.
 * The {@code ActivityProfile} eventually produces an AssemblyRegion based on the {@code AssemblyRegionArgumentCollection} parameters.
 * Let's suppose that this initial region has primary span 350-450 and padded span 100 - 700.
 *
 * Next, the assembly engine assembles all reads that overlap the padded interval to find variant haplotypes and the variants
 * they contain.  The AssemblyRegion is then trimmed down to a new primary interval bound by all assembled variants within the original primary interval
 * and a new padded interval.  The amount of padding of the new padded interval around the variants depends on the needs of local realignment
 * and as such need not equal the original padding that was used for assembly.
 *
 * Re-implementation of the GATK code base. Original author unknown
 * Rust implementation:
 * @author Rhys Newell <rhys.newell@hdr.qut.edu.au>
 */
#[derive(Debug, Clone)]
pub struct AssemblyRegion {
    ref_idx: usize,
    pub(crate) tid: usize,
    contig_length: usize,
    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    pub(crate) reads: Vec<BirdToolRead>,
    /**
     * The active span in which this AssemblyRegion is responsible for calling variants
     */
    pub(crate) active_span: SimpleInterval,
    /**
     * The padded span in which we perform assembly etc in order to call variants within the active span
     */
    pub(crate) padded_span: SimpleInterval,
    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    is_active: bool,
    /**
     * Indicates whether the region has been finalized
     */
    has_been_finalized: bool,
    /**
     * The number of bases within the region that showed large evidence of activity
     */
    activity_density: f32,
}

impl AssemblyRegion {
    /**
     * Create a new AssemblyRegion containing no reads
     *  @param activeSpan the span of this active region
     * @param isActive indicates whether this is an active region, or an inactive one
     * @param padding the active region padding to use for this active region
     */
    pub fn new(
        active_span: SimpleInterval,
        is_active: bool,
        padding: usize,
        contig_length: usize,
        tid: usize,
        ref_idx: usize,
        activity_density: f32
    ) -> AssemblyRegion {
        AssemblyRegion {
            padded_span: AssemblyRegion::make_padded_span(
                &active_span,
                padding,
                contig_length,
                tid,
            ),
            active_span,
            is_active,
            contig_length,
            tid,
            ref_idx,
            reads: Vec::new(),
            has_been_finalized: false,
            activity_density
        }
    }

    pub fn calculate_coverage<L: Locatable>(&self, reads: &[L]) -> f64 {
        let mut coverage = vec![0; self.padded_span.size()];

        for read in reads.iter() {
            let read_start = read.get_start() - self.padded_span.start;
            let read_end = read.get_end() - self.padded_span.start;

            for i in read_start..read_end {
                coverage[i] += 1;
            }
        }

        // calculate the mean
        let total = coverage.iter().sum::<u32>();
        let mean = total as f64 / coverage.len() as f64;
        mean
    }

    pub fn compute_additional_kmer_sizes(&self, current_kmer_sizes: &[usize]) -> Option<Vec<usize>> {
        if self.activity_density < MINIMUM_ACTIVITY_DENSITY_THRESHOLD {
            return None;
        }

        let mut additional_kmer_sizes = Vec::new();

        if (self.activity_density - MINIMUM_ACTIVITY_DENSITY_THRESHOLD) > 0.4 {
            for kmer_size in DEFAULT_ADDITIONAL_KMERS.iter() {
                Self::add_kmer_size(*kmer_size, current_kmer_sizes, &mut additional_kmer_sizes);
            }
        } else if (self.activity_density - MINIMUM_ACTIVITY_DENSITY_THRESHOLD) > 0.2 {
            for kmer_size in DEFAULT_ADDITIONAL_KMERS[0..2].iter() {
                Self::add_kmer_size(*kmer_size, current_kmer_sizes, &mut additional_kmer_sizes);
            }
        } else {
            Self::add_kmer_size(DEFAULT_ADDITIONAL_KMERS[1], current_kmer_sizes, &mut additional_kmer_sizes);
        }

        Some(additional_kmer_sizes)
    }

    fn add_kmer_size(mut kmer_size: usize, current_kmer_sizes: &[usize], additional_kmer_sizes: &mut Vec<usize>) {
        while current_kmer_sizes.contains(&kmer_size) {
            kmer_size += 3;
        }

        additional_kmer_sizes.push(kmer_size);
    }

    pub fn clone_without_reads(&self) -> AssemblyRegion {
        AssemblyRegion {
            padded_span: self.padded_span.clone(),
            active_span: self.active_span.clone(),
            is_active: self.is_active,
            contig_length: self.contig_length,
            tid: self.tid,
            ref_idx: self.ref_idx,
            reads: Vec::new(),
            has_been_finalized: self.has_been_finalized,
            activity_density: self.activity_density
        }
    }

    fn make_padded_span(
        active_span: &SimpleInterval,
        padding: usize,
        contig_length: usize,
        tid: usize,
    ) -> SimpleInterval {
        let mut start = active_span.get_start();
        if padding > start {
            start = 0
        } else {
            start = start - padding
        }

        IntervalUtils::trim_interval_to_contig(
            tid,
            start,
            active_span.get_end() + padding,
            contig_length,
        )
        .unwrap_or_else(|| active_span.clone())
    }

    pub fn new_with_padded_span(
        active_span: SimpleInterval,
        padded_span: SimpleInterval,
        is_active: bool,
        contig_length: usize,
        tid: usize,
        ref_idx: usize,
        activity_density: f32
    ) -> AssemblyRegion {
        AssemblyRegion {
            padded_span,
            active_span,
            is_active,
            contig_length,
            tid,
            ref_idx,
            reads: Vec::new(),
            has_been_finalized: false,
            activity_density
        }
    }

    pub fn get_contig(&self) -> usize {
        self.tid
    }

    pub fn get_start(&self) -> usize {
        self.active_span.get_start()
    }

    pub fn get_end(&self) -> usize {
        self.active_span.get_end()
    }

    pub fn is_active(&self) -> bool {
        self.is_active
    }

    pub fn len(&self) -> usize {
        self.reads.len()
    }

    /**
     * Override activity state of the region
     *
     * Note: Changing the isActive state after construction is a debug-level operation that only engine classes
     * like AssemblyRegionWalker should be able to do
     *
     * @param value new activity state of this region
     */
    // fn set_is_active(&mut self, value: bool) {
    //     self.is_active = value
    // }

    /**
     * Get the span of this assembly region including the padding value
     * @return a non-null SimpleInterval
     */
    pub fn get_padded_span(&self) -> SimpleInterval {
        self.padded_span.clone()
    }

    /**
     * Get the raw span of this assembly region (excluding the padding)
     * @return a non-null SimpleInterval
     */
    pub fn get_span(&self) -> &SimpleInterval {
        &self.active_span
    }

    /**
     * Get an immutable reference of the list of reads currently in this assembly region.
     *
     * The reads are sorted by their coordinate position.
     * @return an unmodifiable and inmutable copy of the reads in the assembly region.
     */
    pub fn get_reads(&self) -> &Vec<BirdToolRead> {
        &self.reads
    }

    /**
     * Move the reads out of an assembly region without losing other information
     * Returns a vector of BirdToolReads
     */
    pub fn move_reads(&mut self) -> Vec<BirdToolRead> {
        let mut container = Vec::with_capacity(self.reads.len());
        std::mem::swap(&mut self.reads, &mut container); // Indiana jones maneuver
        return container;
    }

    /**
     * Get an mutable reference of the list of reads currently in this assembly region.
     *
     * The reads are sorted by their coordinate position.
     * @return an unmodifiable and inmutable copy of the reads in the assembly region.
     */
    pub fn get_reads_mut(&mut self) -> &mut Vec<BirdToolRead> {
        &mut self.reads
    }

    /**
     * Get a clone of the list of reads currently in this assembly region.
     *
     * The reads are sorted by their coordinate position.
     * @return an unmodifiable and inmutable copy of the reads in the assembly region.
     */
    pub fn get_reads_cloned(&self) -> Vec<BirdToolRead> {
        self.reads.clone()
    }

    pub fn take_reads(&mut self) -> Vec<BirdToolRead> {
        std::mem::replace(&mut self.reads, Vec::new())
    }

    /**
     * Trim this region to just the span, producing a new assembly region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param padding the padding size we want for the newly trimmed active region
     * @return a non-null, empty assembly region
     */
    pub fn trim_with_padding(self, span: SimpleInterval, padding: usize) -> AssemblyRegion {
        let padded_span = span.expand_within_contig(padding, self.contig_length);
        return self.trim_with_padded_span(span, padded_span);
    }

    /**
     * Trim this region to no more than the span, producing a new assembly region with properly trimmed reads that
     * attempts to provide the best possible representation of this region covering the span.
     *
     * The challenge here is that span may (1) be larger than can be represented by this assembly region
     * + its original padding and (2) the padding must be symmetric on both sides.  This algorithm
     * therefore determines how best to represent span as a subset of the span of this
     * region with a padding value that captures as much of the span as possible.
     *
     * For example, suppose this active region is
     *
     * Active:    100-200 with padding of 50, so that the true span is 50-250
     * NewExtent: 150-225 saying that we'd ideally like to just have bases 150-225
     *
     * Here we represent the assembly region as a region from 150-200 with 25 bp of padding.
     *
     * The overall constraint is that the region can never exceed the original region, and
     * the padding is chosen to maximize overlap with the desired region
     *
     * @param span the new extend of the active region we want
     * @return a non-null, empty active region
     */
    pub fn trim_with_padded_span(
        self,
        span: SimpleInterval,
        padded_span: SimpleInterval,
    ) -> AssemblyRegion {
        let new_active_span = self.get_span().intersect(&span);
        let new_padded_span = self.get_padded_span().intersect(&padded_span);
        // debug!("new padded span {:?}", &new_padded_span);
        let mut result = AssemblyRegion::new_with_padded_span(
            new_active_span,
            new_padded_span.clone(),
            self.is_active,
            self.contig_length,
            self.tid,
            self.ref_idx,
            self.activity_density
        );

        let mut trimmed_reads = self
            .reads
            .into_iter()
            .map(|read| {
                ReadClipper::hard_clip_to_region(
                    read,
                    new_padded_span.get_start(),
                    new_padded_span.get_end(),
                )
            })
            .filter(|read| !read.is_empty() && read.overlaps(&result.padded_span))
            .collect::<Vec<BirdToolRead>>();
        trimmed_reads.par_sort_unstable();

        result.reads.clear();
        result.reads.extend(trimmed_reads);

        return result;
    }

    /**
     * Add read to this region
     *
     * Read must have alignment start >= than the last read currently in this active region.
     *
     * @throws IllegalArgumentException if read doesn't overlap the padded region of this active region
     *
     * @param read a non-null GATKRead
     */
    pub fn add(&mut self, read: BirdToolRead) {
        // let read_loc = SimpleInterval::new(read.get_contig(), read.get_start(), read.get_end());

        assert!(
            self.padded_span.overlaps(&read),
            "Read does not overlap with active region padded span"
        );

        if !self.reads.is_empty() {
            let final_read = &self.reads.last().unwrap();

            assert!(
                final_read.get_contig() == read.get_contig(),
                "Attempting to add a read to ActiveRegion not on the same contig as other reads"
            );
            assert!(
                read.get_start() >= final_read.get_start(),
                "Attempting to add a read to ActiveRegion out of order w.r.t. other reads: previous {:?} -> current {:?}",
                final_read, &read
            );
        }

        self.reads.push(read);
    }

    pub fn is_valid_read(
        final_read: Option<&BirdToolRead>,
        new_read: &BirdToolRead,
        contig: usize,
    ) -> bool {
        match final_read {
            None => new_read.get_contig() == contig,
            Some(final_read) => {
                return if final_read.get_contig() == new_read.get_contig() {
                    if new_read.get_start() >= final_read.get_start() {
                        true
                    } else {
                        false
                    }
                } else {
                    false
                };
            }
        }
    }

    /**
     * Clear all of the reads currently in this region
     */
    pub fn clear_reads(&mut self) {
        self.reads.clear();
    }

    /**
     * Remove all of the reads in readsToRemove from this region
     * @param readsToRemove the set of reads we want to remove
     */
    pub fn remove_all(mut self, reads_to_remove: &Vec<BirdToolRead>) -> AssemblyRegion {
        self.reads = self
            .reads
            .into_iter()
            .filter(|read| !reads_to_remove.contains(read))
            .collect::<Vec<BirdToolRead>>();

        return self;
    }

    pub fn remove_all_by_reference(
        mut self,
        reads_to_remove: Vec<&BirdToolRead>,
    ) -> AssemblyRegion {
        self.reads = self
            .reads
            .into_iter()
            .filter(|read| !reads_to_remove.contains(&read))
            .collect::<Vec<BirdToolRead>>();

        return self;
    }

    /**
     * Add all readsToAdd to this region
     * @param readsToAdd a collection of readsToAdd to add to this active region
     */
    pub fn add_all(&mut self, reads_to_add: Vec<BirdToolRead>) {
        self.par_extend(reads_to_add)
    }

    /**
     * Get the reference bases from referenceReader spanned by the padded span of this region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases.
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region padded span
     * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    fn get_reference<'a>(
        &self,
        reference_reader: &'a mut ReferenceReader,
        padding: usize,
        genome_loc: &SimpleInterval,
        sequence_already_read_in: bool,
    ) -> &'a [u8] {
        assert!(
            genome_loc.size() > 0,
            "genome_loc must have size > 0 but got {:?}",
            genome_loc
        );

        let interval = SimpleInterval::new(
            self.tid,
            genome_loc.get_start().saturating_sub(padding),
            std::cmp::min(
                *reference_reader
                    .target_lens
                    .get(&self.tid)
                    .unwrap_or(&std::u64::MAX) as usize
                    - 1,
                genome_loc.get_end() + padding,
            ),
        );

        if !sequence_already_read_in {
            reference_reader.update_current_sequence_without_capacity();
            // Update all contig information

            // debug!("Fetching interval... {:?}", &interval);
            reference_reader.fetch_reference_context(self.ref_idx, &interval);
            reference_reader.read_sequence_to_vec();
            reference_reader.current_sequence =
                reference_reader.current_sequence.to_ascii_uppercase();
            if reference_reader.current_sequence.is_empty() {
                panic!(
                    "Retrieved sequence appears to be empty ref_idx {} tid {}",
                    self.ref_idx, self.tid
                );
            };

            return &reference_reader.current_sequence;
        } else {
            return &reference_reader.current_sequence[interval.start
                ..min(
                    interval.get_end() + 1,
                    *reference_reader
                        .target_lens
                        .get(&interval.get_contig())
                        .unwrap() as usize,
                )];
        }
    }

    /**
     * Get the reference bases from referenceReader spanned by the padded span of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region padded region
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    pub fn get_assembly_region_reference<'a>(
        &self,
        reference_reader: &'a mut ReferenceReader,
        padding: usize,
        sequence_already_read_in: bool,
    ) -> &'a [u8] {
        return self.get_reference(
            reference_reader,
            padding,
            &self.padded_span,
            sequence_already_read_in,
        );
    }

    pub fn set_finalized(&mut self, value: bool) {
        self.has_been_finalized = value
    }

    pub fn is_finalized(&self) -> bool {
        self.has_been_finalized
    }
}

impl ParallelExtend<BirdToolRead> for AssemblyRegion {
    fn par_extend<I>(&mut self, par_iter: I)
    where
        I: IntoParallelIterator<Item = BirdToolRead>,
    {
        let par_iter = par_iter
            .into_par_iter()
            .filter(|read| Self::is_valid_read(self.reads.last(), &read, self.tid))
            .collect::<Vec<BirdToolRead>>();
        self.reads.par_extend(par_iter.into_par_iter());
    }
}
