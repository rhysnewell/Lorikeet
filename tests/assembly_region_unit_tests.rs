#![allow(
    non_upper_case_globals,
    non_snake_case
)]

use lorikeet_genome::assembly::assembly_region::AssemblyRegion;
use lorikeet_genome::processing::lorikeet_engine::ReadType;
use lorikeet_genome::reads::bird_tool_reads::BirdToolRead;
use lorikeet_genome::reference::reference_reader::ReferenceReader;
use lorikeet_genome::reference::reference_reader_utils::{ReferenceReaderUtils, GenomesAndContigs};
use lorikeet_genome::utils::artificial_read_utils::ArtificialReadUtils;
use lorikeet_genome::utils::interval_utils::IntervalUtils;
use lorikeet_genome::utils::simple_interval::{Locatable, SimpleInterval};
use std::cmp::{max, min};

#[macro_use]
extern crate lazy_static;

lazy_static! {
    pub static ref CONTIG_LEN: usize = ReferenceReaderUtils::generate_faidx(b37_reference_20_21)
        .index
        .sequences()[0]
        .len as usize;
}
static b37_reference_20_21: &str = "tests/resources/large/human_g1k_v37.20.21.fasta";

#[test]
fn test_constructor() {
    let loc = SimpleInterval::new(1, 10, 20);
    let ar = AssemblyRegion::new(loc.clone(), true, 2, *CONTIG_LEN, 0, 0);
    assert_eq!(ar.is_active(), true);
    assert_eq!(ar.get_span(), &loc);
    assert_eq!(ar.get_contig(), 0)
}

fn test_creating_assembly_regions(
    loc: SimpleInterval,
    is_active: bool,
    extension: usize,
    reader: &mut ReferenceReader,
) {
    let mut region = AssemblyRegion::new(loc.clone(), is_active, extension, *CONTIG_LEN, 0, 0);
    println!(
        "loc {:?} active {} extension {} Contig Length {}",
        &loc, is_active, extension, *CONTIG_LEN
    );
    let dummy_reads = Vec::new();
    assert!(!region.is_finalized());
    assert_eq!(region.get_span(), &loc);
    assert_eq!(
        region.get_padded_span().get_start(),
        loc.get_start().checked_sub(extension).unwrap_or(0)
    );
    assert_eq!(
        region.get_padded_span().get_end(),
        min(loc.get_end() + extension, *CONTIG_LEN)
    );
    assert_eq!(region.is_active(), is_active);
    assert_eq!(region.get_reads(), &dummy_reads);
    assert_eq!(region.len(), 0);

    println!("First");
    assert_good_reference_getter(
        region
            .get_assembly_region_reference(reader, 0, false)
            .to_vec(),
        region.get_padded_span(),
        0,
        reader,
    );
    println!("Second");
    assert_good_reference_getter(
        region
            .get_assembly_region_reference(reader, 0, false)
            .to_vec(),
        region.get_padded_span(),
        0,
        reader,
    );
    println!("Third");
    assert_good_reference_getter(
        region
            .get_assembly_region_reference(reader, 10, false)
            .to_vec(),
        region.get_padded_span(),
        10,
        reader,
    );

    region.set_finalized(false);
    assert!(!region.is_finalized());
    region.set_finalized(true);
    assert!(region.is_finalized());
    region.set_finalized(false);
    assert!(!region.is_finalized());
}

fn assert_good_reference_getter(
    actual_bytes: Vec<u8>,
    span: SimpleInterval,
    padding: usize,
    reader: &mut ReferenceReader,
) {
    let expected_start = span.get_start().checked_sub(padding).unwrap_or(0);
    let expected_stop = min(span.get_end() + padding, *CONTIG_LEN);
    let expected_span = SimpleInterval::new(0, expected_start, expected_stop);
    reader.fetch_reference_context(0, &expected_span);
    println!("Reading sequence to vec");
    reader.read_sequence_to_vec();
    let expected_bytes = &reader.current_sequence;
    assert_eq!(&actual_bytes, expected_bytes);
}

#[test]
fn make_polling_data() {
    println!("Polling data reading in ref...");
    let mut ref_reader = ReferenceReader::new_with_target_names(
        &Some(b37_reference_20_21.to_string()),
        GenomesAndContigs::new(),
        vec![b"20", b"21"],
    );

    ref_reader.target_lens.insert(0, *CONTIG_LEN as u64);

    for start in vec![0, 10, 100, *CONTIG_LEN - 10, *CONTIG_LEN - 1] {
        for size in vec![1, 10, 100, 1000] {
            for ext in vec![0, 1, 10, 100] {
                for is_active in vec![true, false] {
                    let loc = IntervalUtils::trim_interval_to_contig(
                        0,
                        start,
                        start + size - 1,
                        *CONTIG_LEN,
                    )
                    .unwrap();
                    println!("Testing loc {:?}", &loc);
                    test_creating_assembly_regions(loc, is_active, ext, &mut ref_reader)
                }
            }
        }
    }
}

fn test_assembly_region_reads(loc: SimpleInterval, read: BirdToolRead) {
    let mut region = AssemblyRegion::new(loc.clone(), true, 0, *CONTIG_LEN, 0, 0);
    let dummy_reads = Vec::new();
    let dummy_reads_with_one = vec![read.clone()];

    assert_eq!(region.get_reads(), &dummy_reads);
    assert_eq!(region.len(), 0);
    assert_eq!(&region.get_padded_span(), &loc);

    region.add(read.clone());
    assert_eq!(region.get_reads(), &dummy_reads_with_one);
    assert_eq!(region.len(), 1);
    assert_eq!(&region.get_padded_span(), &loc);

    region.clear_reads();
    assert_eq!(region.get_reads(), &dummy_reads);
    assert_eq!(region.len(), 0);
    assert_eq!(&region.get_padded_span(), &loc);

    region.add_all(vec![read.clone()]);
    assert_eq!(region.get_reads(), &dummy_reads_with_one);
    assert_eq!(region.len(), 1);
    assert_eq!(&region.get_padded_span(), &loc);

    let mut region = region.remove_all(&dummy_reads);
    assert_eq!(region.get_reads(), &dummy_reads_with_one);
    assert_eq!(region.len(), 1);
    assert_eq!(&region.get_padded_span(), &loc);

    let mut region = region.remove_all(&dummy_reads_with_one);
    assert_eq!(region.get_reads(), &dummy_reads);
    assert_eq!(region.len(), 0);
    assert_eq!(&region.get_padded_span(), &loc);

    // let read2 = read.clone();
    // read2.read.s
}

#[test]
fn make_assembly_region_reads() {
    println!("Assembly region reading in ref...");

    let mut ref_reader = ReferenceReader::new_with_target_names(
        &Some(b37_reference_20_21.to_string()),
        GenomesAndContigs::new(),
        vec![b"20", b"21"],
    );

    for start in vec![0, 10, 100, *CONTIG_LEN as i64 - 10, *CONTIG_LEN as i64 - 1] {
        for read_start_offset in vec![-100, -10, 0, 10, 100] {
            for read_size in vec![10, 100, 1000] {
                let loc = IntervalUtils::trim_interval_to_contig(
                    0,
                    start as usize,
                    (start + 10) as usize,
                    *CONTIG_LEN,
                )
                .unwrap();

                println!("Make assembly region test {:?}", &loc);
                let read_start = max(start + read_start_offset, 0);
                let read_stop = min(read_start + read_size, *CONTIG_LEN as i64);
                let read_length = read_stop - read_start + 1;
                if read_length > 0 {
                    let read = BirdToolRead::new(
                        ArtificialReadUtils::create_artificial_read_default(
                            "read",
                            0,
                            read_start,
                            read_length as usize,
                            false,
                        ),
                        0,
                        ReadType::Short,
                    );
                    let read_loc = SimpleInterval::new(0, read.get_start(), read.get_end());
                    if read_loc.overlaps(&loc) {
                        test_assembly_region_reads(loc, read)
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------------------------
//
// Make sure bad inputs are properly detected
//
// -----------------------------------------------------------------------------------------------

fn test_bad_read(read1: BirdToolRead, read2: BirdToolRead) {
    let loc = SimpleInterval::new(0, read1.get_start(), read1.get_end());
    let mut region = AssemblyRegion::new(loc, true, 0, *CONTIG_LEN, 0, 0);
    region.add(read1);
    region.add(read2);
}

#[test]
#[should_panic]
fn make_bad_reads_test_1() {
    test_bad_read(
        BirdToolRead::new(
            ArtificialReadUtils::create_artificial_read_default("read1", 0, 10, 10, false),
            0,
            ReadType::Short,
        ),
        BirdToolRead::new(
            ArtificialReadUtils::create_artificial_read_default("read2", 0, 9, 10, false),
            0,
            ReadType::Short,
        ),
    )
}

#[test]
#[should_panic]
fn make_bad_reads_test_2() {
    test_bad_read(
        BirdToolRead::new(
            ArtificialReadUtils::create_artificial_read_default("read1", 0, 10, 10, false),
            0,
            ReadType::Short,
        ),
        BirdToolRead::new(
            ArtificialReadUtils::create_artificial_read_default("read2", 1, 9, 10, false),
            0,
            ReadType::Short,
        ),
    )
}

#[test]
#[should_panic]
fn make_bad_reads_test_3() {
    test_bad_read(
        BirdToolRead::new(
            ArtificialReadUtils::create_artificial_read_default("read1", 1, 10, 10, false),
            0,
            ReadType::Short,
        ),
        BirdToolRead::new(
            ArtificialReadUtils::create_artificial_read_default("read2", 0, 9, 10, false),
            0,
            ReadType::Short,
        ),
    )
}

fn test_trim_assembly_region(
    region_loc: SimpleInterval,
    extension: usize,
    desired_span: SimpleInterval,
    expected_assembly_region: SimpleInterval,
    expected_extension: usize,
) {
    let mut region = AssemblyRegion::new(region_loc, true, extension, *CONTIG_LEN, 0, 0);
    let trimmed = region.trim_with_padded_span(desired_span.clone(), desired_span.clone());
    assert_eq!(
        trimmed.get_span(),
        &expected_assembly_region,
        "Incorrect region"
    );
}

#[test]
fn make_trim_assembly_region_data() {
    // fully enclosed within active region
    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 15, 16),
        SimpleInterval::new(0, 15, 16),
        0,
    );

    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 10, 15),
        SimpleInterval::new(0, 10, 15),
        0,
    );

    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 15, 20),
        SimpleInterval::new(0, 15, 20),
        0,
    );

    // needs extra padding on the right
    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 15, 25),
        SimpleInterval::new(0, 15, 20),
        5,
    );

    // needs extra padding on the left
    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 5, 15),
        SimpleInterval::new(0, 10, 15),
        5,
    );

    // needs extra padding on both
    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 7, 21),
        SimpleInterval::new(0, 10, 20),
        3,
    );
    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 9, 23),
        SimpleInterval::new(0, 10, 20),
        3,
    );

    // desired span captures everything, so we're returning everything.  Tests that extension is set correctly
    test_trim_assembly_region(
        SimpleInterval::new(0, 10, 20),
        10,
        SimpleInterval::new(0, 1, 50),
        SimpleInterval::new(0, 10, 20),
        10,
    );

    // At the start of the chromosome, potentially a bit weird
    test_trim_assembly_region(
        SimpleInterval::new(0, 1, 10),
        10,
        SimpleInterval::new(0, 1, 50),
        SimpleInterval::new(0, 1, 10),
        10,
    );
}
