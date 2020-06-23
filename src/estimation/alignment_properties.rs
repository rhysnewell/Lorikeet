// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

// Taken from Varlociraptor, All credit goes to Johannes and David. Just needed to change
// Bam reader types and generate alignment properties per read to avoid reader being used up

use std::cmp;
use std::error::Error;
use std::f64;
use std::io;
use std::str::FromStr;
use std::u32;

use csv;
use itertools::Itertools;
use ordered_float::NotNan;
use coverm::bam_generator::*;
use rust_htslib::bam::{self, record::Cigar};
use statrs::statistics::{OrderStatistics, Statistics};
use rayon::prelude::*;

use crate::model::variants::Variant;

#[derive(Clone, Debug, Copy, Deserialize, Serialize)]
pub struct AlignmentProperties {
    insert_size: InsertSize,
    pub(crate) max_del_cigar_len: u32,
    pub(crate) max_ins_cigar_len: u32,
    pub(crate) frac_max_softclip: f64,
}

impl AlignmentProperties {
    /// Constructs a dummy instance where all bases are feasible.
    pub fn default(insert_size: InsertSize) -> Self {
        AlignmentProperties {
            insert_size,
            max_del_cigar_len: 30,
            max_ins_cigar_len: 30,
            frac_max_softclip: 1.0,
        }
    }

    /// Update maximum observed cigar operation lengths. Return whether any D, I, S, or H operation
    /// was found in the cigar string.
    pub fn update_max_cigar_ops_len(&mut self, record: &bam::Record) -> bool {
        let norm = |j| NotNan::new(j as f64 / record.seq().len() as f64).unwrap();

        let mut is_regular = true;
        for c in record.cigar().iter() {
            match c {
                &Cigar::SoftClip(l) => {
                    let s = norm(l);
                    self.frac_max_softclip =
                        *cmp::max(s, NotNan::new(self.frac_max_softclip).unwrap());
                    is_regular = false;
                }
                &Cigar::Del(l) => {
                    self.max_del_cigar_len = cmp::max(l, self.max_del_cigar_len);
                    is_regular = false;
                }
                &Cigar::Ins(l) => {
                    self.max_ins_cigar_len = cmp::max(l, self.max_ins_cigar_len);
                    is_regular = false;
                }
                &Cigar::HardClip(_) => {
                    is_regular = false;
                }
                _ => continue,
            }
        }

        is_regular
    }

    /// Estimate `AlignmentProperties` from first 10000 fragments of bam file.
    /// Only reads that are mapped, not duplicates and where quality checks passed are taken.
    pub fn estimate<R: NamedBamReader>(bam: &mut R) -> Result<Self, Box<dyn Error>> {
        let mut properties = AlignmentProperties {
            insert_size: InsertSize::default(),
            max_del_cigar_len: 0,
            max_ins_cigar_len: 0,
            frac_max_softclip: 0.0,
        };

        let mut record = bam::Record::new();
        let mut tlens = Vec::new();
        let mut max_mapq = 0;
        let mut i = 0;
        while i <= 10000 {
            if !bam.read(&mut record)? {
                break;
            }

            if record.is_unmapped() || record.is_duplicate() || record.is_quality_check_failed() {
                continue;
            }

            max_mapq = cmp::max(max_mapq, record.mapq());

            let is_regular = properties.update_max_cigar_ops_len(&record);

            if is_regular
                && !record.is_mate_unmapped()
                && record.is_first_in_template()
                && record.tid() == record.mtid()
                && record.mapq() > 0
            {
                // record insert size
                tlens.push(record.insert_size().abs() as f64);
            }

            i += 1;
        }

        let upper = tlens.percentile(95);
        let lower = tlens.percentile(5);
        let mut valid = tlens
            .into_iter()
            .filter(|l| *l <= upper && *l >= lower)
            .collect_vec();
        properties.insert_size.mean = valid.median();
        properties.insert_size.sd = valid.iter().std_dev();

        Ok(properties)
    }

    /// Estimate `AlignmentProperties` from each record in bam file.
    /// Only reads that are mapped, not duplicates and where quality checks passed are taken.
    /// Return will need to be extended to existing vector
    pub fn estimate_from_record(record: &bam::Record, properties: &mut AlignmentProperties) -> Option<f64> {

        let mut tlen = None;
        if record.is_unmapped() || record.is_duplicate() || record.is_quality_check_failed() {
            tlen
        } else {
            let is_regular = properties.update_max_cigar_ops_len(&record);

            if is_regular
                && !record.is_mate_unmapped()
                && record.is_first_in_template()
                && record.tid() == record.mtid()
                && record.mapq() > 0
            {
                // record insert size
                tlen = Some(record.insert_size().abs() as f64);
            }
            tlen
        }
    }

    pub fn update_properties(tlens: &mut Vec<f64>, properties: &mut AlignmentProperties) -> Result<Self, Box<dyn Error>> {
        let upper = tlens.percentile(95);
        let lower = tlens.percentile(5);
        let mut valid: Vec<_> = tlens
            .into_par_iter()
            .filter(|l| {
                **l <= upper && **l >= lower
            }).map(|l| *l)
            .collect();
        properties.insert_size.mean = valid.median();
        properties.insert_size.sd = valid.iter().std_dev();
        Ok(*properties)
    }

    /// Number of bases that are feasible for overlapping the variant.
    pub fn feasible_bases(&self, read_len: u32, variant: &Variant) -> u32 {
        match variant {
            &Variant::Deletion(l) if l <= self.max_del_cigar_len => read_len,
            &Variant::Insertion(ref seq) if seq.len() as u32 <= self.max_ins_cigar_len => read_len,
            &Variant::SNV(_) => return read_len,
            &Variant::None => return read_len,
            _ => (read_len as f64 * self.frac_max_softclip) as u32,
        }
    }

    pub fn insert_size(&self) -> &InsertSize {
        &self.insert_size
    }
}

/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub struct InsertSize {
    pub mean: f64,
    pub sd: f64,
}

impl InsertSize {
    /// Obtain insert size from samtools stats output.
    pub fn from_samtools_stats<R: io::Read>(
        samtools_stats: &mut R,
    ) -> Result<InsertSize, Box<dyn Error>> {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .flexible(true)
            .from_reader(samtools_stats);

        let mut insert_size = InsertSize::default();

        for rec in rdr.records() {
            let rec = rec?;
            if &rec[0] == "SN" {
                if &rec[1] == "insert size average:" {
                    insert_size.mean = f64::from_str(&rec[2])?;
                } else if &rec[1] == "insert size standard deviation:" {
                    insert_size.sd = f64::from_str(&rec[2])?;
                    break;
                }
            }
        }

        Ok(insert_size)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use coverm::*;
    use rust_htslib::bam::record::Record;

    #[test]
    fn test_estimate() {
        let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
            vec!["tests/data/test2.bam"]);
        for reader in bam_readers {
            let mut bam = reader.start();

            let props = AlignmentProperties::estimate(&mut bam).unwrap();
            println!("{:?}", props);

            assert_eq!(props.insert_size.mean, 499.0);
            assert_eq!(props.insert_size.sd, 39.541883346487914);
            assert_eq!(props.max_del_cigar_len, 6);
            assert_eq!(props.max_ins_cigar_len, 5);
            assert_eq!(props.frac_max_softclip, 0.44666666666666666);
        }
    }

    #[test]
    fn test_estimate_from_records() {
        let bam_readers = bam_generator::generate_named_bam_readers_from_bam_files(
            vec!["tests/data/test2.bam"]);
        for reader in bam_readers {
            let mut bam = reader.start();
            let mut record= Record::new();
            let mut tlens = Vec::new();
            let mut properties = AlignmentProperties {
                insert_size: InsertSize::default(),
                max_del_cigar_len: 0,
                max_ins_cigar_len: 0,
                frac_max_softclip: 0.0,
            };

            while bam.read(&mut record)
                .expect("Error while reading BAM record") == true {
                match AlignmentProperties::estimate_from_record(&mut record, &mut properties) {
                    Some(tlen) => tlens.push(tlen),
                    None => {},
                };
            }

            let props = AlignmentProperties::update_properties(&mut tlens, &mut properties).unwrap();
            println!("{:?}", props);

            assert_eq!(props.insert_size.mean, 499.0);
            assert_eq!(props.insert_size.sd, 39.541883346487914);
            assert_eq!(props.max_del_cigar_len, 6);
            assert_eq!(props.max_ins_cigar_len, 5);
            assert_eq!(props.frac_max_softclip, 0.44666666666666666);
        }
    }

    #[test]
    fn test_parse_insert_size() {
        let insert_size = InsertSize::from_samtools_stats(&mut io::BufReader::new(
            fs::File::open("tests/resources/samtools_stats.example.txt").unwrap(),
        ))
        .unwrap();
        assert_eq!(insert_size.mean, 311.7);
        assert_eq!(insert_size.sd, 15.5);
    }
}
