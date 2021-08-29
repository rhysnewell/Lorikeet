use std::io::Error;

#[derive(Debug, Clone)]
pub enum BirdToolError {
    InvalidClip(String),
    CigarBuilderError(String),
    InvalidLocation(String),
    NonContiguousIntervals(String),
}
