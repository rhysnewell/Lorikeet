use std::error::Error;
use std::fmt;

#[derive(Clone)]
pub enum BirdToolError {
    InvalidClip(String),
    CigarBuilderError(String),
    IOError(String),
    InvalidLocation(String),
    NonContiguousIntervals(String),
    SkipException(String),
    InvalidVariationEvent(String),
    ProcessPanicked(String),
    DebugError(String),
}

// Implement std::fmt::Display for AppError
impl fmt::Display for BirdToolError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "An Error Occurred, Please Try Again!") // user-facing output
    }
}

// Implement std::fmt::Debug for AppError
impl fmt::Debug for BirdToolError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{ file: {}, line: {} }}", file!(), line!()) // programmer-facing output
    }
}

impl Error for BirdToolError {
    fn description(&self) -> &str {
        match self {
            BirdToolError::InvalidClip(val)
            | BirdToolError::IOError(val)
            | BirdToolError::CigarBuilderError(val)
            | BirdToolError::InvalidLocation(val)
            | BirdToolError::NonContiguousIntervals(val)
            | BirdToolError::SkipException(val)
            | BirdToolError::InvalidVariationEvent(val)
            | BirdToolError::ProcessPanicked(val)
            | BirdToolError::DebugError(val) => val,
        }
    }
}
