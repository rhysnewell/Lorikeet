pub mod bam_generator;
pub mod mapping_index_maintenance;
pub mod mapping_parameters;
pub mod filter;

use rust_htslib::bam::record::Record;
use std::sync::Arc;

use crate::parse_percentage;


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
        debug!("Generating output writer from path {:?}", file_path_opt);
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
                debug!("Outputing to STDOUT");
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

fn aux_as(record: &rust_htslib::bam::Record) -> i64 {
    match record.aux("AS".as_bytes()) {
        Ok(value) => {
            if let rust_htslib::bam::record::Aux::U8(v) = value {
                v as i64
            } else if let rust_htslib::bam::record::Aux::U16(v) = value {
                v as i64
            } else {
                panic!("Unexpected data type of AS aux tag, found {:?}", value)
            }
        }
        Err(e) => {
            panic!(
                "Mapping record encountered that does not have an 'AS' \
                        auxiliary tag in the SAM/BAM format. This is required \
                        for ranking pairs of alignments. Error was {}",
                e
            )
        }
    }
}

#[derive(Debug)]
struct FilterParameters {
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
}
impl FilterParameters {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> FilterParameters {
        let f = FilterParameters {
            flag_filters: FlagFilter {
                include_improper_pairs: !m.get_flag("proper-pairs-only"),
                include_secondary: m.get_flag("include-secondary"),
                include_supplementary: !m.get_flag("exclude-supplementary"),
            },
            min_aligned_length_single: *m.get_one::<u32>("min-read-aligned-length").unwrap_or(&0),
            min_percent_identity_single: parse_percentage(&m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(&m, "min-read-aligned-percent"),
            min_aligned_length_pair: *m
                .get_one::<u32>("min-read-aligned-length-pair")
                .unwrap_or(&0),
            min_percent_identity_pair: parse_percentage(&m, "min-read-percent-identity-pair"),
            min_aligned_percent_pair: parse_percentage(&m, "min-read-aligned-percent-pair"),
        };
        debug!("Filter parameters set as {:?}", f);
        return f;
    }

    // pub fn add_metabat_filtering_if_required(&mut self, _m: &clap::ArgMatches) {
    //     // if doing_metabat(&m) {
    //     //     info!(
    //     //         "Setting single read percent identity threshold at 0.97 for \
    //     //          MetaBAT adjusted coverage, and not filtering out supplementary, \
    //     //          secondary and improper pair alignments"
    //     //     );
    //     //     // we use >= where metabat uses >. Gah.
    //     //     self.min_percent_identity_single = 0.97001;
    //     //     self.flag_filters.include_improper_pairs = true;
    //     //     self.flag_filters.include_supplementary = true;
    //     //     self.flag_filters.include_secondary = true;
    //     // }
    // }

    pub fn doing_filtering(&self) -> bool {
        return self.min_percent_identity_single > 0.0
            || self.min_percent_identity_pair > 0.0
            || self.min_aligned_percent_single > 0.0
            || self.min_aligned_percent_pair > 0.0
            || self.min_aligned_length_single > 0
            || self.min_aligned_length_pair > 0;
    }
}
