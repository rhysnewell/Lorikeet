use model::allele_frequency_calculator::AlleleFrequencyCalculator;

pub struct GenotypingEngine {
    allele_frequency_calculator: AlleleFrequencyCalculator,
    number_of_genomes: usize,
    samples: Vec<String>,
    do_allele_specific_calcs: bool,
}

impl GenotypingEngine {
    pub fn new(
        samples: Vec<String>,
        do_allele_specific_calcs: bool,
    )
}