pub struct VCFConstants {}

impl VCFConstants {
    pub const GENOTYPE_KEY: String = "GT".to_string();
    pub const GENOTYPE_POSTERIORS_KEY: String = "GP".to_string();
    pub const GENOTYPE_PRIOR_KEY: String = "PG".to_string();
    pub const LOW_QUAL_FILTER_NAME: String = "LowQual".to_string();
}