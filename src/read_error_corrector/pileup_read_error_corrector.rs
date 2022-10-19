use read_error_corrector::read_error_corrector::ReadErrorCorrector;
use reads::bird_tool_reads::BirdToolRead;

pub struct PileupReadErrorCorrector {
    log_odds_threshold: f64,
}

impl PileupReadErrorCorrector {
    pub fn new(log_odds_threshold: f64) -> Self {
        Self { log_odds_threshold }
    }
}
//
// impl ReadErrorCorrector for PileupReadErrorCorrector {
//     fn correct_reads(&mut self, reads: Vec<BirdToolRead>) -> Vec<BirdToolRead> {
//
//     }
// }
