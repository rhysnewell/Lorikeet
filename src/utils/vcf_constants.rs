lazy_static! {
    pub static ref GENOTYPE_KEY: String = "GT".to_string();
    pub static ref GENOTYPE_POSTERIORS_KEY: String = "GP".to_string();
    pub static ref GENOTYPE_PRIOR_KEY: String = "PG".to_string();
    pub static ref LOW_QUAL_FILTER_NAME: String = "LowQual".to_string();
    pub static ref MLE_ALLELE_COUNT_KEY: String = "MLEAC".to_string();
    pub static ref MLE_ALLELE_FREQUENCY_KEY: String = "MLEAF".to_string();
    pub static ref AS_QUAL_KEY: String = "AS_QUAL".to_string();
    pub static ref NUMBER_OF_DISCOVERED_ALLELES_KEY: String = "NDA".to_string();
}

pub struct VCFConstants {}

impl VCFConstants {
    pub const SPANNING_DELETION_ALLELE: char = '*';
    pub const NO_CALL_ALLELE: char = '.';
    pub const NULL_ALLELE: char = '-';
}
