use std::io::Read;
use std::sync::atomic::{compiler_fence, Ordering};
use nix::sys::stat;
use nix::unistd;
use rust_htslib::errors::Result as HtslibResult;
use std::path::Path;
use tempdir::TempDir;
use rust_htslib::bam;
use rust_htslib::bam::{FetchDefinition, Read as BamRead};
use rust_htslib::errors::Error;

use crate::bam_parsing::FlagFilter;
use crate::bam_parsing::filter::ReferenceSortedBamFilter;
use crate::bam_parsing::mapping_index_maintenance::MappingIndex;
use crate::bam_parsing::mapping_parameters::ReadFormat;

use tempfile;

pub trait NamedBamReader {
    // Name of the stoit
    fn name(&self) -> &str;

    // Read a record into record parameter
    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>>;

    // Return pileup alignments
    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::Reader>>;

    // Return the bam header of the final BAM file
    fn header(&self) -> &bam::HeaderView;

    fn path(&self) -> &str;

    fn finish(self);

    //set the number of threads for Bam reading
    fn set_threads(&mut self, n_threads: usize);

    // Number of reads that were detected
    fn num_detected_primary_alignments(&self) -> u64;
}

pub trait NamedBamReaderGenerator<T> {
    // For readers that map, start the process of mapping
    fn start(self) -> T;
}

pub trait IndexedNamedBamReader {
    // Name of the stoit
    fn name(&self) -> &str;

    // Fetch the specified region
    fn fetch<'a, T: Into<FetchDefinition<'a>>>(&mut self, fetch_definition: T)
        -> Result<(), Error>;

    // Read a record into record parameter
    fn read(&mut self, record: &mut bam::record::Record) -> bool;

    // Return pileup alignments
    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::IndexedReader>>;

    // Return the bam header of the final BAM file
    fn header(&self) -> &bam::HeaderView;

    fn path(&self) -> &str;

    fn finish(self);

    //set the number of threads for Bam reading
    fn set_threads(&mut self, n_threads: usize);

    // Number of reads that were detected
    fn num_detected_primary_alignments(&self) -> u64;
}

#[derive(Debug, Clone, Copy)]
#[allow(non_camel_case_types)]
pub enum MappingProgram {
    BWA_MEM,
    BWA_MEM2,
    MINIMAP2_SR,
    MINIMAP2_ONT,
    MINIMAP2_HIFI,
    MINIMAP2_PB,
    MINIMAP2_NO_PRESET,
}

pub struct BamFileNamedReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    num_detected_primary_alignments: u64,
    path: String,
}

#[derive(Debug)]
pub struct IndexedBamFileNamedReader {
    stoit_name: String,
    bam_reader: bam::IndexedReader,
    num_detected_primary_alignments: u64,
    path: String,
}

impl NamedBamReader for BamFileNamedReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }

    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        let res = self.bam_reader.read(record);
        if res == Some(Ok(())) && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        return res;
    }

    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::Reader>> {
        Some(self.bam_reader.pileup())
    }

    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }

    fn path(&self) -> &str {
        &self.path
    }

    fn finish(self) {}

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.bam_reader.set_threads(n_threads - 1).unwrap();
        }
    }

    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments;
    }
}

impl NamedBamReaderGenerator<BamFileNamedReader> for BamFileNamedReader {
    fn start(self) -> BamFileNamedReader {
        BamFileNamedReader {
            stoit_name: self.stoit_name,
            bam_reader: self.bam_reader,
            num_detected_primary_alignments: 0,
            path: self.path,
        }
    }
}

impl IndexedNamedBamReader for IndexedBamFileNamedReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }

    // Fetch the specified region
    fn fetch<'a, T: Into<FetchDefinition<'a>>>(
        &mut self,
        fetch_definition: T,
    ) -> Result<(), Error> {
        self.bam_reader.fetch(fetch_definition)
    }

    fn read(&mut self, record: &mut bam::record::Record) -> bool {
        let res = self.bam_reader.read(record);
        let res = match res {
            Some(result) => match result {
                Ok(_) => {
                    if !record.is_secondary() && !record.is_supplementary() {
                        self.num_detected_primary_alignments += 1;
                    };
                    true
                }
                Err(e) => panic!("Error: {:?}", e),
            },
            None => false,
        };
        return res;
    }

    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::IndexedReader>> {
        Some(self.bam_reader.pileup())
    }

    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }

    fn path(&self) -> &str {
        &self.path
    }

    fn finish(self) {}

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.bam_reader.set_threads(n_threads - 1).unwrap();
        }
    }

    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments;
    }
}

impl NamedBamReaderGenerator<IndexedBamFileNamedReader> for IndexedBamFileNamedReader {
    fn start(self) -> IndexedBamFileNamedReader {
        IndexedBamFileNamedReader {
            stoit_name: self.stoit_name,
            bam_reader: self.bam_reader,
            num_detected_primary_alignments: 0,
            path: self.path,
        }
    }
}

pub struct StreamingNamedBamReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    tempdir: TempDir,
    path: String,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    num_detected_primary_alignments: u64,
    minimap2_log_file_index: Option<usize>,
}

pub struct StreamingNamedBamReaderGenerator {
    stoit_name: String,
    tempdir: TempDir,
    cache_path: String,
    fifo_path: std::path::PathBuf,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    minimap2_log_file_index: Option<usize>,
}

impl NamedBamReaderGenerator<StreamingNamedBamReader> for StreamingNamedBamReaderGenerator {
    fn start(self) -> StreamingNamedBamReader {
        // debug!("Starting mapping processes");
        let mut processes = vec![];
        // let mut i = 0;
        for mut preprocess in self.pre_processes {
            // debug!("Running mapping command: {}", self.command_strings[i]);
            // i += 1;
            processes.push(preprocess.spawn().expect("Unable to execute bash"));
        }
        let bam_reader = match bam::Reader::from_path(&self.fifo_path) {
            Ok(reader) => reader,
            Err(upstream_error) => {
                error!(
                    "Failed to correctly find or parse BAM file at {:?}: {}",
                    self.fifo_path, upstream_error
                );
                complete_processes(
                    processes,
                    self.command_strings,
                    self.log_file_descriptions,
                    self.log_files,
                    Some(self.tempdir),
                );
                panic!("Failure to find or parse BAM file, cannot continue");
            }
        };
        return StreamingNamedBamReader {
            stoit_name: self.stoit_name,
            bam_reader: bam_reader,
            tempdir: self.tempdir,
            path: self.cache_path,
            processes: processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
            num_detected_primary_alignments: 0,
            minimap2_log_file_index: self.minimap2_log_file_index,
        };
    }
}

pub fn complete_processes(
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    _tempdir: Option<TempDir>,
) {
    let mut failed_any = false;
    let mut overall_stderrs = vec![];
    for mut process in processes {
        // debug!("Finishing process...");
        let es = process
            .wait()
            .expect("Failed to glean exitstatus from mapping process");
        // debug!("Finshed {:?}", es.success());
        let failed = !es.success();
        if failed || log_enabled!(log::Level::Debug) {
            if failed {
                failed_any = true;
                error!("Error when running mapping process. Exitstatus was {:?}. Command run was: {:?}", es, command_strings);
            } else {
                // debug!("Successfully finished process {:?}", process);
            }
            let mut err = String::new();
            process
                .stderr
                .expect("Failed to grab stderr from failed mapping process")
                .read_to_string(&mut err)
                .expect("Failed to read stderr into string");
            // debug!("The overall STDERR was: {:?}", err);
            overall_stderrs.push(err);
        }
    }
    if failed_any || log_enabled!(log::Level::Debug) {
        for (description, tf) in log_file_descriptions.iter().zip(log_files.into_iter()) {
            let mut contents = String::new();
            tf.into_file()
                .read_to_string(&mut contents)
                .expect(&format!("Failed to read log file for {}", description));
            if failed_any {
                error!("The STDERR for the {:} part was: {}", description, contents);
            } else {
                // debug!("The STDERR for the {:} part was: {}", description, contents);
            }
        }
    }
    if failed_any {
        panic!("Cannot continue since mapping failed.");
    }
    // debug!("Process finished without error.");

    // There's (maybe) a (difficult to reproduce) single-thread issue where the
    // tempdir gets dropped before the process is finished. Hopefully putting a
    // compiler fence here stops this.
    compiler_fence(Ordering::SeqCst);
    // debug!("After fence, for tempdir {:?}", tempdir);
}

impl NamedBamReader for StreamingNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }

    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        let res = self.bam_reader.read(record);
        if res == Some(Ok(())) && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        return res;
    }

    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::Reader>> {
        Some(self.bam_reader.pileup())
    }

    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }

    fn path(&self) -> &str {
        &self.path
    }

    fn finish(self) {
        // Check minimap2 didn't complain about unequal numbers of reads
        match self.minimap2_log_file_index {
            None => {}
            Some(log_file_index) => {
                let mut contents = String::new();
                std::fs::File::open(&self.log_files[log_file_index])
                    .expect("Failed to read minimap2 log file")
                    .read_to_string(&mut contents)
                    .expect("Failed to read minimap2 log file to string");
                if contents.contains("query files have different number of records") {
                    error!("The STDERR for the minimap2 part was: {}", contents);
                    panic!(
                        "Not continuing since when input file pairs have \
                        unequal numbers of reads this usually means \
                        incorrect / corrupt files were specified"
                    );
                }
            }
        };

        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            Some(self.tempdir),
        );
    }

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.bam_reader.set_threads(n_threads - 1).unwrap();
        }
    }

    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments;
    }
}

pub fn generate_named_bam_readers_from_bam_files(bam_paths: Vec<&str>) -> Vec<BamFileNamedReader> {
    bam_paths
        .iter()
        .map(|path| BamFileNamedReader {
            stoit_name: std::path::Path::new(path)
                .file_stem()
                .unwrap()
                .to_str()
                .expect("failure to convert bam file name to stoit name - UTF8 error maybe?")
                .to_string(),
            bam_reader: bam::Reader::from_path(path)
                .expect(&format!("Unable to find BAM file {}", path)),
            num_detected_primary_alignments: 0,
            path: path.to_string(),
        })
        .collect()
}

pub fn generate_indexed_named_bam_readers_from_bam_files(
    bam_paths: Vec<&str>,
    threads: u32,
) -> Vec<IndexedBamFileNamedReader> {
    bam_paths
        .iter()
        .map(|path| {
            // check and build bam index if it doesn't exist
            if !Path::new(&(path.to_string() + ".bai")).exists() {
                bam::index::build(
                    path,
                    Some(&format!("{}.bai", path).as_str()),
                    bam::index::Type::Bai,
                    threads,
                )
                .expect(&format!("Unable to index bam at {}", &path));
            }
            IndexedBamFileNamedReader {
                stoit_name: std::path::Path::new(path)
                    .file_stem()
                    .unwrap()
                    .to_str()
                    .expect("failure to convert bam file name to stoit name - UTF8 error maybe?")
                    .to_string(),
                bam_reader: bam::IndexedReader::from_path(path)
                    .expect(&format!("Unable to find BAM file {}", path)),
                num_detected_primary_alignments: 0,
                path: path.to_string(),
            }
        })
        .collect()
}

pub fn generate_named_bam_readers_from_reads(
    mapping_program: MappingProgram,
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    discard_unmapped: bool,
    mapping_options: Option<&str>,
    include_reference_in_stoit_name: bool,
) -> StreamingNamedBamReaderGenerator {
    let tmp_dir = TempDir::new("coverm_fifo").expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let mapping_log = tempfile::NamedTempFile::new().expect(&format!(
        "Failed to create {:?} log tempfile",
        mapping_program
    ));
    let samtools2_log =
        tempfile::NamedTempFile::new().expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log =
        tempfile::NamedTempFile::new().expect("Failed to create cache samtools view log tempfile");

    let cached_bam_file_args = match cached_bam_file {
        Some(path) => {
            format!(
                "|tee {} |samtools view {} -@ {} -b -o '{}' 2>{}",
                // tee
                fifo_path.as_os_str().to_str().unwrap(),
                // samtools view
                match discard_unmapped {
                    true => "-F4",
                    false => "",
                },
                threads - 1,
                path,
                samtools_view_cache_log
                    .path()
                    .to_str()
                    .expect("Failed to convert tempfile path to str")
            )
        }
        None => format!("> {}", fifo_path.as_os_str().to_str().unwrap()),
    };

    let mapping_command = build_mapping_command(
        mapping_program,
        read_format,
        threads,
        read1_path,
        reference,
        read2_path,
        mapping_options,
    );
    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile_in(tmp_dir.path())
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         {} 2>{} \
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         {}",
        // Mapping program
        mapping_command,
        mapping_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix
            .path()
            .to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads - 1,
        samtools2_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // Caching (or not)
        cached_bam_file_args
    );
    // debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    // Required because of https://github.com/wwood/CoverM/issues/58
    let minimap2_log_file_index = match mapping_program {
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => None,
        // Required because of https://github.com/lh3/minimap2/issues/527
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_NO_PRESET => Some(0),
    };

    let mut log_descriptions = vec![
        format!("{:?}", mapping_program).to_string(),
        "samtools sort".to_string(),
    ];
    let mut log_files = vec![mapping_log, samtools2_log];
    if cached_bam_file.is_some() {
        log_descriptions.push("samtools view for cache".to_string());
        log_files.push(samtools_view_cache_log);
    }

    let stoit_name = match include_reference_in_stoit_name {
        true => {
            std::path::Path::new(reference)
                .file_name()
                .expect("Unable to convert reference to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string()
                + "/"
        }
        false => "".to_string(),
    } + &std::path::Path::new(read1_path)
        .file_name()
        .expect("Unable to convert read1 name to file name")
        .to_str()
        .expect("Unable to covert file name into str")
        .to_string();

    return StreamingNamedBamReaderGenerator {
        stoit_name: stoit_name,
        tempdir: tmp_dir,
        fifo_path: fifo_path,
        cache_path: match cached_bam_file {
            Some(cache_path) => cache_path.to_string(),
            None => "".to_string(),
        },
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c \"{}\"", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files: log_files,
        minimap2_log_file_index: minimap2_log_file_index,
    };
}

pub struct FilteredBamReader {
    stoit_name: String,
    filtered_stream: ReferenceSortedBamFilter,
    path: String,
}

impl NamedBamReader for FilteredBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }

    fn read(&mut self, mut record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        self.filtered_stream.read(&mut record)
    }

    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::Reader>> {
        Some(self.filtered_stream.pileup())
    }

    fn header(&self) -> &bam::HeaderView {
        &self.filtered_stream.reader.header()
    }

    fn path(&self) -> &str {
        &self.path
    }

    fn finish(self) {}
    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.filtered_stream
                .reader
                .set_threads(n_threads - 1)
                .unwrap();
        }
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        return self.filtered_stream.num_detected_primary_alignments;
    }
}

impl NamedBamReaderGenerator<FilteredBamReader> for FilteredBamReader {
    fn start(self) -> FilteredBamReader {
        FilteredBamReader {
            stoit_name: self.stoit_name,
            filtered_stream: self.filtered_stream,
            path: self.path,
        }
    }
}

pub fn generate_filtered_bam_readers_from_bam_files(
    bam_paths: Vec<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
) -> Vec<FilteredBamReader> {
    let mut generators: Vec<FilteredBamReader> = vec![];

    for path in bam_paths {
        let filtered: FilteredBamReader;
        let stoit_name = std::path::Path::new(path)
            .file_stem()
            .unwrap()
            .to_str()
            .expect("failure to convert bam file name to stoit name - UTF8 error maybe?")
            .to_string();
        let reader =
            bam::Reader::from_path(path).expect(&format!("Unable to find BAM file {}", path));

        filtered = FilteredBamReader {
            stoit_name: stoit_name,
            filtered_stream: ReferenceSortedBamFilter::new(
                reader,
                flag_filters.clone(),
                min_aligned_length_single,
                min_percent_identity_single,
                min_aligned_percent_single,
                min_aligned_length_pair,
                min_percent_identity_pair,
                min_aligned_percent_pair,
                true,
            ),
            path: path.to_string(),
        };

        generators.push(filtered)
    }

    return generators;
}

pub struct StreamingFilteredNamedBamReader {
    stoit_name: String,
    filtered_stream: ReferenceSortedBamFilter,
    tempdir: TempDir,
    path: String,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

pub struct StreamingFilteredNamedBamReaderGenerator {
    pub stoit_name: String,
    pub tempdir: TempDir,
    pub fifo_path: std::path::PathBuf,
    pub cache_path: String,
    pub pre_processes: Vec<std::process::Command>,
    pub command_strings: Vec<String>,
    pub flag_filters: FlagFilter,
    pub min_aligned_length_single: u32,
    pub min_percent_identity_single: f32,
    pub min_aligned_percent_single: f32,
    pub min_aligned_length_pair: u32,
    pub min_percent_identity_pair: f32,
    pub min_aligned_percent_pair: f32,
    pub log_file_descriptions: Vec<String>,
    pub log_files: Vec<tempfile::NamedTempFile>,
}

impl NamedBamReaderGenerator<StreamingFilteredNamedBamReader>
    for StreamingFilteredNamedBamReaderGenerator
{
    fn start(self) -> StreamingFilteredNamedBamReader {
        // debug!("Starting mapping processes");
        let mut processes = vec![];
        for mut preprocess in self.pre_processes {
            processes.push(preprocess.spawn().expect("Unable to execute bash"));
        }
        let bam_reader = match bam::Reader::from_path(&self.fifo_path) {
            Ok(reader) => reader,
            Err(upstream_error) => {
                error!(
                    "Failed to correctly find or parse BAM file at {:?}: {}",
                    self.fifo_path, upstream_error
                );
                complete_processes(
                    processes,
                    self.command_strings,
                    self.log_file_descriptions,
                    self.log_files,
                    Some(self.tempdir),
                );
                panic!("Failure to find or parse BAM file, cannot continue");
            }
        };

        let filtered_stream = ReferenceSortedBamFilter::new(
            bam_reader,
            self.flag_filters,
            self.min_aligned_length_single,
            self.min_percent_identity_single,
            self.min_aligned_percent_single,
            self.min_aligned_length_pair,
            self.min_percent_identity_pair,
            self.min_aligned_percent_pair,
            true,
        );
        return StreamingFilteredNamedBamReader {
            stoit_name: self.stoit_name,
            filtered_stream: filtered_stream,
            tempdir: self.tempdir,
            path: self.cache_path,
            processes: processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
        };
    }
}

impl NamedBamReader for StreamingFilteredNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }

    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        self.filtered_stream.read(record)
    }

    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::Reader>> {
        Some(self.filtered_stream.pileup())
    }

    fn header(&self) -> &bam::HeaderView {
        self.filtered_stream.reader.header()
    }

    fn path(&self) -> &str {
        &self.path
    }

    fn finish(self) {
        // debug!(
        //     "Finishing StreamingFilteredNamedBamReader. Tempdir is {:?}",
        //     self.tempdir.path()
        // );
        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            Some(self.tempdir),
        )
    }

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.filtered_stream
                .reader
                .set_threads(n_threads - 1)
                .unwrap();
        }
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        return self.filtered_stream.num_detected_primary_alignments;
    }
}

pub fn generate_filtered_named_bam_readers_from_reads(
    mapping_program: MappingProgram,
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    bwa_options: Option<&str>,
    discard_unmapped: bool,
    include_reference_in_stoit_name: bool,
) -> StreamingFilteredNamedBamReaderGenerator {
    let streaming = generate_named_bam_readers_from_reads(
        mapping_program,
        reference,
        read1_path,
        read2_path,
        read_format,
        threads,
        cached_bam_file,
        discard_unmapped,
        bwa_options,
        include_reference_in_stoit_name,
    );
    return StreamingFilteredNamedBamReaderGenerator {
        stoit_name: streaming.stoit_name,
        tempdir: streaming.tempdir,
        fifo_path: streaming.fifo_path,
        cache_path: streaming.cache_path,
        pre_processes: streaming.pre_processes,
        command_strings: streaming.command_strings,
        log_file_descriptions: streaming.log_file_descriptions,
        log_files: streaming.log_files,
        flag_filters: flag_filters,
        min_aligned_length_single: min_aligned_length_single,
        min_percent_identity_single: min_percent_identity_single,
        min_aligned_percent_single: min_aligned_percent_single,
        min_aligned_length_pair: min_aligned_length_pair,
        min_percent_identity_pair: min_percent_identity_pair,
        min_aligned_percent_pair: min_aligned_percent_pair,
    };
}

pub struct BamGeneratorSet<T> {
    pub generators: Vec<T>,
    pub index: Option<Box<dyn MappingIndex>>,
}

pub struct NamedBamMaker {
    stoit_name: String,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

pub struct NamedBamMakerGenerator {
    stoit_name: String,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

pub fn generate_bam_maker_generator_from_reads(
    mapping_program: MappingProgram,
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: &str,
    discard_unmapped: bool,
    mapping_options: Option<&str>,
) -> NamedBamMakerGenerator {
    let mapping_log = tempfile::NamedTempFile::new().expect(&format!(
        "Failed to create {:?} log tempfile",
        mapping_program
    ));
    let samtools2_log =
        tempfile::NamedTempFile::new().expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log =
        tempfile::NamedTempFile::new().expect("Failed to create cache samtools view log tempfile");

    let mapping_command = build_mapping_command(
        mapping_program,
        read_format,
        threads,
        read1_path,
        reference,
        read2_path,
        mapping_options,
    );

    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile()
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         {} 2>{} \
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         | samtools view {} -b -@ {} -o '{}' 2>{}",
        // Mapping program
        mapping_command,
        mapping_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix
            .path()
            .to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads - 1,
        samtools2_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // samtools view
        match discard_unmapped {
            true => "-F4",
            false => "",
        },
        threads - 1,
        cached_bam_file,
        samtools_view_cache_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str")
    );
    // debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    let log_descriptions = vec![
        format!("{:?}", mapping_program).to_string(),
        "samtools sort".to_string(),
        "samtools view for cache".to_string(),
    ];
    let log_files = vec![mapping_log, samtools2_log, samtools_view_cache_log];

    return NamedBamMakerGenerator {
        stoit_name: std::path::Path::new(reference)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string()
            + "/"
            + &std::path::Path::new(read1_path)
                .file_name()
                .expect("Unable to convert read1 name to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string(),
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c \"{}\"", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files: log_files,
    };
}

impl NamedBamReaderGenerator<NamedBamMaker> for NamedBamMakerGenerator {
    fn start(self) -> NamedBamMaker {
        // debug!("Starting mapping processes");
        let mut processes = vec![];
        // let mut i = 0;
        for mut preprocess in self.pre_processes {
            // debug!("Running mapping command: {}", self.command_strings[i]);
            // i += 1;
            processes.push(preprocess.spawn().expect("Unable to execute bash"));
        }
        return NamedBamMaker {
            stoit_name: self.stoit_name,
            processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
        };
    }
}

impl NamedBamMaker {
    pub fn name(&self) -> &str {
        &(self.stoit_name)
    }

    pub fn complete(&self) {}

    pub fn finish(self) {
        // debug!("Finishing NamedBamMaker...");
        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            None,
        )
    }
}

pub fn build_mapping_command(
    mapping_program: MappingProgram,
    read_format: ReadFormat,
    threads: u16,
    read1_path: &str,
    reference: &str,
    read2_path: Option<&str>,
    mapping_options: Option<&str>,
) -> String {
    let read_params1 = match mapping_program {
        // minimap2 auto-detects interleaved based on read names
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => "",
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => match read_format {
            ReadFormat::Interleaved => "-p",
            ReadFormat::Coupled | ReadFormat::Single => "",
        },
    };

    let read_params2 = match read_format {
        ReadFormat::Interleaved => format!("'{}'", read1_path),
        ReadFormat::Coupled => format!("'{}' '{}'", read1_path, read2_path.unwrap()),
        ReadFormat::Single => format!("'{}'", read1_path),
    };

    return format!(
        "{} {} -t {} {} '{}' {}",
        match mapping_program {
            MappingProgram::BWA_MEM => "bwa mem".to_string(),
            MappingProgram::BWA_MEM2 => "bwa-mem2 mem".to_string(),
            _ => {
                let split_prefix = tempfile::Builder::new()
                    .prefix("coverm-minimap2-split-index")
                    .tempfile()
                    .expect(&format!(
                        "Failed to create {:?} minimap2 split_prefix file",
                        mapping_program
                    ));
                format!(
                    "minimap2 --split-prefix {} -a {}",
                    split_prefix
                        .path()
                        .to_str()
                        .expect("Failed to convert split prefix tempfile path to str"),
                    match mapping_program {
                        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => unreachable!(),
                        MappingProgram::MINIMAP2_SR => "-x sr",
                        MappingProgram::MINIMAP2_ONT => "-x map-ont",
                        MappingProgram::MINIMAP2_HIFI => "-x map-hifi",
                        MappingProgram::MINIMAP2_PB => "-x map-pb",
                        MappingProgram::MINIMAP2_NO_PRESET => "",
                    }
                )
            }
        },
        mapping_options.unwrap_or(""),
        threads,
        read_params1,
        reference,
        read_params2
    );
}

pub struct PlaceholderBamFileReader {
    header: bam::HeaderView,
}

impl NamedBamReader for PlaceholderBamFileReader {
    fn name(&self) -> &str {
        &("placeholder")
    }

    fn read(&mut self, _record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        None
    }

    fn pileup(&mut self) -> Option<bam::pileup::Pileups<bam::Reader>> {
        None
    }

    fn header(&self) -> &bam::HeaderView {
        &self.header
    }

    fn path(&self) -> &str {
        &("placeholder")
    }

    fn finish(self) {}

    fn set_threads(&mut self, _n_threads: usize) {}

    fn num_detected_primary_alignments(&self) -> u64 {
        0
    }
}

impl NamedBamReaderGenerator<PlaceholderBamFileReader> for PlaceholderBamFileReader {
    fn start(self) -> PlaceholderBamFileReader {
        PlaceholderBamFileReader {
            header: bam::HeaderView::from_header(&bam::Header::new()),
        }
    }
}

pub fn generate_placeholder() -> Vec<PlaceholderBamFileReader> {
    vec![PlaceholderBamFileReader {
        header: bam::HeaderView::from_header(&bam::Header::new()),
    }]
}
