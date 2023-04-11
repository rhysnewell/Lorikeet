pub mod bam_generator;
pub mod mapping_index_maintenance;
pub mod mapping_parameters;
pub mod filter;

use rust_htslib::bam::record::Record;
use std::sync::Arc;


pub const CONCATENATED_FASTA_FILE_SEPARATOR: &str = "~";

#[derive(PartialEq, Debug, Clone)]
pub struct ReadsMapped {
    num_mapped_reads: u64,
    num_reads: u64,
}

#[derive(Clone, Debug)]
pub struct FlagFilter {
    pub include_improper_pairs: bool,
    pub include_supplementary: bool,
    pub include_secondary: bool,
}

impl FlagFilter {
    pub fn passes(&self, record: &Record) -> bool {
        if !self.include_secondary && record.is_secondary() {
            return false;
        }
        if !self.include_supplementary && record.is_supplementary() {
            return false;
        }
        if !self.include_improper_pairs && !record.is_proper_pair() {
            return false;
        }
        return true;
    }
}

pub struct OutputWriter {
    pub output_file: Option<Arc<std::sync::Mutex<std::fs::File>>>,
}

impl OutputWriter {
    pub fn generate(file_path_opt: Option<&str>) -> OutputWriter {
        // debug!("Generating output writer from path {:?}", file_path_opt);
        match file_path_opt {
            Some(file_path) => {
                if file_path == "-" {
                    info!("Outputing to STDOUT");
                    OutputWriter { output_file: None }
                } else {
                    let path = std::path::Path::new(&file_path);
                    info!("Writing output to file: {}", file_path);
                    OutputWriter {
                        output_file: Some(Arc::new(std::sync::Mutex::new(
                            std::fs::File::create(path)
                                .expect(&format!("Failed to create output file: {}", file_path)),
                        ))),
                    }
                }
            }
            None => {
                // debug!("Outputing to STDOUT");
                OutputWriter { output_file: None }
            }
        }
    }
}

impl std::io::Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match &mut self.output_file {
            None => std::io::stdout().write(buf),
            Some(f) => f.lock().expect("failed to unlock output file").write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match &mut self.output_file {
            None => std::io::stdout().flush(),
            Some(f) => f.lock().expect("failed to unlock output file").flush(),
        }
    }
}

impl Clone for OutputWriter {
    fn clone(&self) -> OutputWriter {
        OutputWriter {
            output_file: match &self.output_file {
                Some(f) => Some(f.clone()),
                None => None,
            },
        }
    }
}

// Verbose to do this many times throughout the code so making a function to
// abstract.
fn nm(record: &rust_htslib::bam::Record) -> u64 {
    match record.aux("NM".as_bytes()) {
        Ok(value) => {
            if let rust_htslib::bam::record::Aux::U8(v) = value {
                v as u64
            } else if let rust_htslib::bam::record::Aux::U16(v) = value {
                v as u64
            } else {
                panic!("Unexpected data type of NM aux tag, found {:?}", value)
            }
        }
        Err(e) => {
            panic!(
                "Mapping record encountered that does not have an 'NM' \
                        auxiliary tag in the SAM/BAM format. This is required \
                        to work out some coverage statistics. Error was {}",
                e
            )
        }
    }
}
