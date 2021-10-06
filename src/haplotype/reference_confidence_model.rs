use assembly::assembly_region::AssemblyRegion;
use haplotype::haplotype::Haplotype;
use model::byte_array_allele::Allele;
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
    pub fn create_reference_haplotype<L: Locatable>(
        active_region: &AssemblyRegion,
        ref_bases: &[u8],
        padded_reference_loc: &SimpleInterval,
    ) -> Haplotype<L> {
        debug!(
            "Active region span {:?} padded ref loc {:?} start {} -> {}",
            active_region.get_padded_span(),
            padded_reference_loc,
            active_region.get_padded_span().get_start(),
            padded_reference_loc.get_start()
        );
        let alignment_start = active_region.get_padded_span().get_start() as i64
            - padded_reference_loc.get_start() as i64;
        if alignment_start < 0 {
            panic!("Bad alignment start for reference haplotype creation");
        }
        let mut ref_haplotype = Haplotype::new(ref_bases, true);
        ref_haplotype.set_alignment_start_hap_wrt_ref(alignment_start as usize);
        let mut c = Vec::new();
        c.push(Cigar::Match(ref_haplotype.get_bases().len() as u32));
        ref_haplotype.set_cigar(c);

        return ref_haplotype;
    }
}
