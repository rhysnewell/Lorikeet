use assembly::assembly_region::AssemblyRegion;
use haplotype::haplotype::Haplotype;
use rust_htslib::bam::record::Cigar;
use utils::simple_interval::{Locatable, SimpleInterval};

pub struct ReferenceConfidenceModel {}

impl ReferenceConfidenceModel {
    /**
     * Create a reference haplotype for an active region
     *
     * @param activeRegion the active region
     * @param refBases the ref bases
     * @param paddedReferenceLoc the location spanning of the refBases -- can be longer than activeRegion.getLocation()
     * @return a reference haplotype
     */
    pub fn create_reference_haplotype<'a, L: Locatable>(
        active_region: &AssemblyRegion,
        ref_bases: &[u8],
        padded_reference_loc: &SimpleInterval,
    ) -> Haplotype<'a, L> {
        let alignment_start = active_region
            .get_padded_span()
            .get_start()
            .checked_sub(padded_reference_loc.get_start())
            .unwrap_or(panic!("Bad alignment start in create_reference_haplotype"));

        let mut ref_haplotype = Haplotype::new(ref_bases, true);
        ref_haplotype.set_alignment_start_hap_wrt_ref(alignment_start);
        let mut c = Vec::new();
        c.push(Cigar::Match(ref_haplotype.get_bases().len() as u32));
        ref_haplotype.set_cigar(c);

        return ref_haplotype;
    }
}
