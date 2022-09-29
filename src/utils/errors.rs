
#[derive(Debug, Clone)]
pub enum BirdToolError {
    InvalidClip(String),
    CigarBuilderError(String),
    InvalidLocation(String),
    NonContiguousIntervals(String),
    SkipException(String),
    InvalidVariationEvent(String),
    ProcessPanicked(String),
    DebugError(String)
}
