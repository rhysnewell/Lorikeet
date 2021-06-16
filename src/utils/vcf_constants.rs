pub struct VCFConstants {}

impl VCFConstants {
    pub const GENOTYPE_KEY: String = "GT".to_string();
    pub const GENOTYPE_POSTERIORS_KEY: String = "GP".to_string();
    pub const GENOTYPE_PRIOR_KEY: String = "PG".to_string();
    pub const LOW_QUAL_FILTER_NAME: String = "LowQual".to_string();
    pub const MLE_ALLELE_COUNT_KEY: String = "MLEAC".to_string();
    pub const MLE_ALLELE_FREQUENCY_KEY: String = "MLEAF".to_string();
    pub const AS_QUAL_KEY: String = "AS_QUAL".to_string();
    pub const NUMBER_OF_DISCOVERED_ALLELES_KEY: String = "NDA".to_string();
}