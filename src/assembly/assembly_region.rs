use rust_htslib::bam::Record;
use utils::simple_interval::SimpleInterval;
use utils::interval_utils::IntervalUtils;


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
 */
pub struct AssemblyRegion {
    tid: usize,
    contig_length: usize,
    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    reads: Vec<Read>,
    /**
     * The active span in which this AssemblyRegion is responsible for calling variants
     */
    active_span: SimpleInterval,
    /**
     * The padded span in which we perform assembly etc in order to call variants within the active span
     */
    padded_span: SimpleInterval,
    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    is_active: bool,
    /**
     * Indicates whether the region has been finalized
     */
    has_been_finalized: bool,
}

impl AssemblyRegion {
    /**
     * Create a new AssemblyRegion containing no reads
     *  @param activeSpan the span of this active region
     * @param isActive indicates whether this is an active region, or an inactive one
     * @param padding the active region padding to use for this active region
     */
    pub fn new(active_span: SimpleInterval, is_active: bool, padding: usize, contig_length: usize, tid: usize) -> AssemblyRegion {
        AssemblyRegion {
            padded_span: AssemblyRegion::make_padded_span(&active_span, padding, contig_length, tid),
            active_span,
            is_active,
            contig_length,
            tid,
            reads: Vec::new(),
            has_been_finalized: false,
        }
    }

    fn make_padded_span(active_span: &SimpleInterval, padding: usize, contig_length: usize, tid: usize) -> SimpleInterval {
        let mut start = active_span.get_start();
        if padding > start {
            start = 0
        } else {
            start = start - padding
        }

        IntervalUtils::trim_interval_to_contig(tid, start, active_span.get_end() + padding, contig_length);
    }

    pub fn new_with_padded_span(
        active_span: SimpleInterval,
        padded_span: SimpleInterval,
        is_active: bool,
        contig_length: usize,
        tid: usize
    ) -> AssemblyRegion {
        AssemblyRegion {
            padded_span,
            active_span,
            is_active,
            contig_length,
            tid,
            reads: Vec::new(),
            has_been_finalized: false,
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

    /**
     * Override activity state of the region
     *
     * Note: Changing the isActive state after construction is a debug-level operation that only engine classes
     * like AssemblyRegionWalker should be able to do
     *
     * @param value new activity state of this region
     */
    fn set_is_active(&mut self, value: bool) {
        self.is_active = value
    }

    /**
     * Get the span of this assembly region including the padding value
     * @return a non-null SimpleInterval
     */
    pub fn get_padded_span(&self) -> &SimpleInterval {
        &self.padded_span
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
    pub fn get_reads(&self) -> &Vec<Read> {
        &self.reads
    }

    /**
     * Trim this region to just the span, producing a new assembly region without any reads that has only
     * the extent of newExtend intersected with the current extent
     * @param span the new extend of the active region we want
     * @param padding the padding size we want for the newly trimmed active region
     * @return a non-null, empty assembly region
     */
    pub fn trim_with_padding(&self, span: SimpleInterval, padding: usize) -> AssemblyRegion {
        let padded_span = span.expand_within_contig(padding, self.contig_length);
        return self.trim_with_padded_span(span, padded_span)
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
    pub fn trim_with_padded_span(&self, span: SimpleInterval, padded_span: SimpleInterval) -> AssemblyRegion {
        let new_active_span = self.get_span().intersect(&span);
        let new_padded_span = self.get_padded_span().intersect(&padded_span);

        let result = AssemblyRegion::new_with_padded_span(
            new_active_span, new_padded_span, self.is_active, self.contig_length, self.tid
        );

        let trimmed_reads = self.reads.
    }
}