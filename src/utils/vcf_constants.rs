use model::byte_array_allele::ByteArrayAllele;
use model::variants::SPAN_DEL_ALLELE;

pub static DEPTH_KEY: &str = "DP";
pub static GENOTYPE_KEY: &str = "GT";
pub static HAPLOTYPE_CALLER_PHASING_ID_KEY: &str = "PID";
pub static HAPLOTYPE_CALLER_PHASING_GT_KEY: &str = "PGT";
pub static PHASE_SET_KEY: &str = "PS";
pub static PHASE_QUALITY_KEY: &str = "PQ";
pub static GENOTYPE_POSTERIORS_KEY: &str = "GP";
pub static GENOTYPE_PRIOR_KEY: &str = "PG";
pub static LOW_QUAL_FILTER_NAME: &str = "LowQual";
pub static MLE_ALLELE_COUNT_KEY: &str = "MLEAC";
pub static MLE_ALLELE_FREQUENCY_KEY: &str = "MLEAF";
pub static AS_QUAL_KEY: &str = "AS_QUAL";
pub static NUMBER_OF_DISCOVERED_ALLELES_KEY: &str = "NDA";

// lazy_static! {
//
//     // pub static ref DEPTH_KEY: String = "DP".to_string();
//     // pub static ref GENOTYPE_KEY: String = "GT".to_string();
//     // pub static ref GENOTYPE_POSTERIORS_KEY: String = "GP".to_string();
//     // pub static ref GENOTYPE_PRIOR_KEY: String = "PG".to_string();
//     // pub static ref LOW_QUAL_FILTER_NAME: String = "LowQual".to_string();
//     // pub static ref MLE_ALLELE_COUNT_KEY: String = "MLEAC".to_string();
//     // pub static ref MLE_ALLELE_FREQUENCY_KEY: String = "MLEAF".to_string();
//     // pub static ref AS_QUAL_KEY: String = "AS_QUAL".to_string();
//     // pub static ref NUMBER_OF_DISCOVERED_ALLELES_KEY: String = "NDA".to_string();
//
//     // FORMAT keys
//
//
// }

pub struct VCFConstants {}

impl VCFConstants {
    pub const SPANNING_DELETION_ALLELE: char = '*';
    pub const NO_CALL_ALLELE: char = '.';
    pub const NULL_ALLELE: char = '-';

    pub fn is_spanning_deletion(allele: &ByteArrayAllele) -> bool {
        return allele == &*SPAN_DEL_ALLELE;
    }
}
