use crate::assembly::assembly_region::AssemblyRegion;
use crate::model::variant_context::VariantContext;
use crate::model::variant_context_utils::VariantContextUtils;
use crate::utils::simple_interval::Locatable;

/**
 * Tandem repeat unit composition and counts per allele
 *
 * <p>This annotation tags variants that fall within tandem repeat sets. It also provides the composition of the tandem repeat units and the number of times they are repeated for each allele (including the REF allele).</p>
 *
 * <p>A tandem repeat unit is composed of one or more nucleotides that are repeated multiple times in series. Repetitive sequences are difficult to map to the reference because they are associated with multiple alignment possibilities. Knowing the number of repeat units in a set of tandem repeats tells you the number of different positions the tandem repeat can be placed in. The observation of many tandem repeat units multiplies the number of possible representations that can be made of the region.
 */
pub struct TandemRepeat {}

impl TandemRepeat {
    pub fn get_num_tandem_repeat_units(
        region: &AssemblyRegion,
        vc: &mut VariantContext,
        reference: &[u8],
    ) -> Option<(Vec<usize>, Vec<u8>)> {
        let start_index = vc.loc.get_start() + 1 - region.get_padded_span().get_start(); // +1 to exclude leading match base common to VC's ref and alt alleles
        return VariantContextUtils::get_num_tandem_repeat_units(
            vc,
            &reference[start_index..reference.len()],
        );
    }
}
