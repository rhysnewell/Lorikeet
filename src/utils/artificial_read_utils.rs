use estimation::lorikeet_engine::ReadType;
use reads::bird_tool_reads::BirdToolRead;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{Read, Reader, Record};
use rust_htslib::{
    bam,
    bam::{Header, HeaderView},
};
use std::rc::Rc;

pub static DEFAULT_READ_GROUP_PREFIX: &str = "ReadGroup";
pub static DEFAULT_PLATFORM_UNIT_PREFIX: &str = "Lane";
pub static DEFAULT_PLATFORM_PREFIX: &str = "Platform";
pub static DEFAULT_SAMPLE_NAME: &str = "SampleX";
pub static DEFAULT_PROGRAM_NAME: &str = "Program";
pub static READ_GROUP_ID: &str = "x";

pub struct ArtificialReadUtils {}

impl ArtificialReadUtils {
    pub const DEFAULT_READ_LENGTH: usize = 50;

    pub fn get_default_reader() -> Reader {
        return bam::Reader::from_path("tests/resources/small/test_finalize_region_reads.sam")
            .unwrap();
    }

    pub fn create_artificial_read(bases: &[u8], qual: &[u8], cigar: CigarString) -> BirdToolRead {
        let bam = Self::get_default_reader();
        // let header = Rc::new(bam.header().clone());
        let mut record = Record::new();

        record.set("default_read".as_bytes(), Some(&cigar), bases, qual);
        record.set_header(Rc::new(bam.header().clone()));
        record.set_pos(10000);
        record.set_tid(0);

        return BirdToolRead::new(record, 0, ReadType::Short);
    }
}
