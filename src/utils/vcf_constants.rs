use crate::model::byte_array_allele::ByteArrayAllele;
use crate::model::variants::SPAN_DEL_ALLELE;

lazy_static! {

    pub static ref DEPTH_KEY: String = "DP".to_string();
    pub static ref GENOTYPE_KEY: String = "GT".to_string();
    pub static ref GENOTYPE_POSTERIORS_KEY: String = "GP".to_string();
    pub static ref GENOTYPE_PRIOR_KEY: String = "PG".to_string();
    pub static ref LOW_QUAL_FILTER_NAME: String = "LowQual".to_string();
    pub static ref MLE_ALLELE_COUNT_KEY: String = "MLEAC".to_string();
    pub static ref MLE_ALLELE_FREQUENCY_KEY: String = "MLEAF".to_string();
    pub static ref AS_QUAL_KEY: String = "AS_QUAL".to_string();
    pub static ref NUMBER_OF_DISCOVERED_ALLELES_KEY: String = "NDA".to_string();
    pub static ref HAPLOTYPE_CALLER_PHASING_ID_KEY: String = "PID".to_string();
    pub static ref HAPLOTYPE_CALLER_PHASING_GT_KEY: String = "PGT".to_string();
    pub static ref PHASE_SET_KEY: String = "PS".to_string();
    pub static ref PHASE_QUALITY_KEY: String = "PQ".to_string();

    // FORMAT keys
    pub static ref STRAND_COUNT_BY_SAMPLE_KEY: String = "SAC".to_string();
    pub static ref PHRED_SCALED_POSTERIORS_KEY: String = "PP".to_string();
}

pub struct VCFConstants {}

impl VCFConstants {
    pub const SPANNING_DELETION_ALLELE: char = '*';
    pub const NO_CALL_ALLELE: char = '.';
    pub const NULL_ALLELE: char = '-';

    pub fn is_spanning_deletion(allele: &ByteArrayAllele) -> bool {
        return allele == &*SPAN_DEL_ALLELE;
    }
}
