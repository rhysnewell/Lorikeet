lazy_static! {
    pub static ref FLAT_BETA: BetaDistributionShape = BetaDistributionShape::new(1.0, 1.0);
}

pub struct BetaDistributionShape {
    alpha: f64,
    beta: f64,
}

impl BetaDistributionShape {
    pub fn new(alpha: f64, beta: f64) -> BetaDistributionShape {
        assert!(
            alpha > 0.0,
            "Alpha must be greater than 0 but got {}",
            alpha
        );
        assert!(beta > 0.0, "Beta must be greater than 0 but got {}", beta);

        BetaDistributionShape { alpha, beta }
    }

    pub fn get_alpha(&self) -> f64 {
        self.alpha
    }

    pub fn get_beta(&self) -> f64 {
        self.beta
    }
}
