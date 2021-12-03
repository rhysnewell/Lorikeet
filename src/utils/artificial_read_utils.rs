use processing::lorikeet_engine::ReadType;
use reads::bird_tool_reads::BirdToolRead;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{Read, Reader, Record};
use rust_htslib::{
    bam,
    bam::{Header, HeaderView},
};
use std::convert::TryFrom;
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

    pub fn create_artificial_read_with_name_and_pos(
        qname: String,
        tid: i32,
        alignment_start: i64,
        bases: &[u8],
        qual: &[u8],
        cigar: &str,
        sample_index: usize,
    ) -> BirdToolRead {
        let bam = Self::get_default_reader();
        // let header = Rc::new(bam.header().clone());
        let mut record = Record::new();

        record.set(
            qname.as_bytes(),
            Some(&CigarString::try_from(cigar).unwrap()),
            bases,
            qual,
        );
        record.set_header(Rc::new(bam.header().clone()));
        record.set_pos(alignment_start);
        record.set_tid(tid);

        return BirdToolRead::new(record, sample_index, ReadType::Short);
    }

    pub fn create_artificial_read_default(
        name: &str,
        tid: usize,
        alignment_start: i64,
        length: usize,
        unmapped: bool,
    ) -> Record {
        let bam = Self::get_default_reader();
        let mut record = Record::new();

        let cigar = CigarString(vec![Cigar::Match(length as u32)]);
        let seq = vec!['A' as u8; length];
        let quals = vec![30; length];

        record.set(
            name.as_bytes(),
            Some(&cigar),
            seq.as_slice(),
            quals.as_slice(),
        );
        record.set_pos(alignment_start);
        record.set_tid(tid as i32);
        record.unset_proper_pair();

        if unmapped {
            record.set_unmapped();
        };

        return record;
    }
}
