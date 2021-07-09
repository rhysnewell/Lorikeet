use reads::bird_tool_reads::BirdToolRead;

pub trait ReadErrorCorrector {
    fn correct_reads(reads: Vec<BirdToolRead>) -> Vec<BirdToolRead>;
}