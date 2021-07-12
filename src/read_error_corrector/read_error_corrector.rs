use reads::bird_tool_reads::BirdToolRead;
use read_error_corrector::nearby_kmer_error_corrector::NearbyKmerErrorCorrector;

pub trait ReadErrorCorrector {
    fn correct_reads(&mut self, reads: &Vec<BirdToolRead>) -> Vec<BirdToolRead>;
}