
use std::fmt::Debug;
use std::ops::Range;
use std::collections::{HashSet, HashMap};

use bio::stats::LogProb;
use itertools::Itertools;
use ordered_float::NotNan;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};
use rust_htslib::{bcf, bcf::record::Numeric};

use rayon::prelude::*;

pub type AlleleFreq = NotNan<f64>;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord)]
pub enum StrandBias {
    None,
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None
    }
}

impl StrandBias {
    pub fn is_some(&self) -> bool {
        if let StrandBias::None = self {
            false
        } else {
            true
        }
    }

    pub fn forward_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_one(),
            StrandBias::Reverse => LogProb::ln_zero(),
        }
    }

    pub fn reverse_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_zero(),
            StrandBias::Reverse => LogProb::ln_one(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, EnumString, EnumIter, IntoStaticStr)]
pub enum VariantType {
    #[strum(serialize = "INS")]
    Insertion(Option<Range<u32>>),
    #[strum(serialize = "DEL")]
    Deletion(Option<Range<u32>>),
    #[strum(serialize = "SNV")]
    SNV,
    #[strum(serialize = "MNV")]
    MNV,
    #[strum(serialize = "REF")]
    None, // site with no suggested alternative allele
}

impl From<&str> for VariantType {
    fn from(string: &str) -> VariantType {
        match string {
            "INS" => VariantType::Insertion(None),
            "DEL" => VariantType::Deletion(None),
            "SNV" => VariantType::SNV,
            "REF" => VariantType::None,
            _ => panic!("bug: given string does not describe a valid variant type"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub enum SVType {
    Dup
}

#[derive(Clone, Copy, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub struct SV {
    sv: SVType,
    len: u32,
    start: u32,
    end: u32,
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub enum Variant {
    Deletion(u32),
    Insertion(Vec<u8>),
    SNV(u8),
    MNV(Vec<u8>),
    SV(SV),
    None,
}

#[allow(unused)]
impl Variant {

    pub fn has_fragment_evidence(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::SV(_) => true,
            &Variant::SNV(_) => false,
            &Variant::MNV(_) => false,
            &Variant::None => false,
        }
    }

    pub fn is_single_base(&self) -> bool {
        match self {
            &Variant::SNV(_) | &Variant::None => true,
            _ => false,
        }
    }

    pub fn is_snv(&self) -> bool {
        match self {
            &Variant::SNV(_) => true,
            _ => false,
        }
    }

    pub fn is_indel(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::SV(_) => false,
            &Variant::SNV(_) => false,
            &Variant::MNV(_) => false,
            &Variant::None => false,
        }
    }

    pub fn is_type(&self, vartype: &VariantType) -> bool {
        match (self, vartype) {
            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
                l >= range.start && l < range.end
            }
            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            }
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::SNV(_), &VariantType::SNV) => true,
            (&Variant::MNV(_), &VariantType::MNV) => true,
            (&Variant::None, &VariantType::None) => true,
            _ => false,
        }
    }

    pub fn end(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length,
            &Variant::Insertion(_) => start + 1, // end of insertion is the next regular base
            &Variant::SV(sv) => sv.end,
            &Variant::SNV(_) | &Variant::None => start,
            &Variant::MNV(ref alt) => start + alt.len() as u32,
        }
    }

    pub fn centerpoint(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length / 2,
            &Variant::Insertion(_) => start, // end of insertion is the next regular base
            &Variant::SV(sv) => (sv.start + sv.len) / 2,
            &Variant::SNV(_) | &Variant::None => start,
            &Variant::MNV(ref alt) => start + alt.len() as u32 / 2,
        }
    }

    pub fn len(&self) -> u32 {
        match self {
            &Variant::Deletion(l) => l,
            &Variant::Insertion(ref s) => s.len() as u32,
            &Variant::SV(sv) => sv.len,
            &Variant::SNV(_) => 1,
            &Variant::MNV(ref alt) => alt.len() as u32,
            &Variant::None => 1,
        }
    }
}

/// The filter tag given to the locus
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum Filter {
    LowCov,
    Amb,
    Del,
    PASS,
    None,
}

#[allow(unused)]
impl Filter {
    pub fn from(string: &str) -> Filter {
        match string {
            "PASS" => Filter::PASS,
            "LowCov" => Filter::LowCov,
            "Amb" => Filter::Amb,
            "Del" => Filter::Del,
            _ => Filter::None,
        }
    }

    pub fn from_result(string: Result<&str, std::str::Utf8Error>) -> Filter {
        match string {
            Ok("PASS") => Filter::PASS,
            Ok("LowCov") => Filter::LowCov,
            Ok("Amb") => Filter::Amb,
            Ok("Del") => Filter::Del,
            _ => Filter::None,
        }
    }
}

/// Information about each base position
#[derive(Clone, Debug, PartialEq)]
pub struct Base {
    // Contig ID
    pub tid: u32,
    // Position on contig 0-based
    pub pos: i64,
    // Reference allele
    pub refr: Vec<u8>,
    // Alternate allele Variant enum
    pub variant: Variant,
    // Filter tag
    pub filters: Vec<HashSet<Filter>>,
    // Depth of good quality estimation.reads
    pub depth: Vec<i32>,
    // Depth including bad quality estimation.reads
    pub truedepth: Vec<i32>,
    // Depth as decided by CoverM
    pub totaldepth: Vec<i32>,
    pub referencedepth: Vec<i32>,
    // Physical coverage of valid inserts across locus
    pub physicalcov: Vec<i32>,
    // Mean base quality at locus
    pub baseq: Vec<i32>,
    // Mean read mapping quality at locus
    pub mapq: Vec<i32>,
    // Variant confidence / quality by depth
    pub conf: Vec<i32>,
    // Nucleotide count at each locus
    pub nucs: HashMap<char, Vec<i32>>,
    // Percentage of As, Cs, Gs, Ts weighted by Q & MQ at locus
    pub pernucs: HashMap<char, Vec<i32>>,
    // insertion count at locus
    pub ic: Vec<i32>,
    // deletion count at locus
    pub dc: Vec<i32>,
    // number of estimation.reads clipped here
    pub xc: Vec<i32>,
    // allele count in genotypes, for each ALT allele.
    pub ac: Vec<i32>,
    // fraction in support for alternate allele
    pub af: Vec<f64>,
    // Frequency of variant
    pub freq: Vec<f64>,
    // Read ids assigned to variant
    pub reads: HashSet<Vec<u8>>,
    // CLR transformed relative abundances
    pub rel_abunds: Vec<f64>
}

#[allow(unused)]
impl Base {
    pub fn add_depth(&mut self, sample_idx: usize, d: i32) {
        if self.totaldepth[sample_idx] == 0 {
            self.totaldepth[sample_idx] = d;
//            if self.variant == Variant::None {
//                self.depth[sample_idx] = d;
//            }
        }
    }

    pub fn combine_sample(&mut self, other: &Base, sample_idx: usize, total_depth: i32) {
        if &self != &other {
            self.filters[sample_idx] = other.filters[sample_idx].clone();
            self.depth[sample_idx] = other.depth[sample_idx];
            self.truedepth[sample_idx] = other.truedepth[sample_idx];
            self.totaldepth[sample_idx] = total_depth;
            self.physicalcov[sample_idx] = other.physicalcov[sample_idx];
            self.baseq[sample_idx] = other.baseq[sample_idx];
            self.mapq[sample_idx] = other.mapq[sample_idx];
            self.conf[sample_idx] = other.conf[sample_idx];
            self.ic[sample_idx] = other.ic[sample_idx];
            self.dc[sample_idx] = other.dc[sample_idx];
            self.xc[sample_idx] = other.xc[sample_idx];
            self.ac[sample_idx] = other.ac[sample_idx];
            self.af[sample_idx] = other.af[sample_idx];
            self.freq[sample_idx] = other.freq[sample_idx];

        } else {
            self.totaldepth[sample_idx] = total_depth;
        }
    }

    pub fn update_total_depth(&mut self, depth: i32, sample_idx: usize) {
        self.totaldepth[sample_idx] = depth
    }

    pub fn new(tid: u32, pos: i64, refr: Vec<u8>, sample_count: usize) -> Base {
        Base {
            tid,
            pos,
            refr,
            variant: Variant::None,
            filters: vec![HashSet::new(); sample_count],
            depth: vec![0; sample_count],
            truedepth: vec![0; sample_count],
            totaldepth: vec![0; sample_count],
            referencedepth: vec![0; sample_count],
            physicalcov: vec![0; sample_count],
            baseq: vec![0; sample_count],
            mapq: vec![0; sample_count],
            conf: vec![0; sample_count],
            nucs: HashMap::new(),
            pernucs: HashMap::new(),
            ic: vec![0; sample_count],
            dc: vec![0; sample_count],
            xc: vec![0; sample_count],
            ac: vec![0; sample_count],
            af: vec![0.; sample_count],
            freq: vec![0.; sample_count],
            rel_abunds: vec![0.; sample_count],
            reads: HashSet::new(),
        }
    }

    pub fn from_vcf_record(record: &mut bcf::Record, sample_count: usize, sample_idx: usize, longread: bool) -> Option<Vec<Base>> {

        let variants = collect_variants(record, false,
                                        false, None);
        if variants.len() > 0 {

            // Separate filters into hashset of filter struct
            let mut filter_hash = HashSet::new();
            {
                let filters = record.filters();
                let header = record.header();
                for filter in filters {
                    filter_hash.insert(Filter::from_result(std::str::from_utf8(&header.id_to_name(filter)[..])));
                }
            }
            let mut bases = vec!();
            let mut refr_base_empty = true;
            for (idx, variant) in variants.iter().enumerate() {
                // Get elements from record
                let mut base = Base::new(record.rid().unwrap(),
                                         record.pos(), record.alleles()[idx].to_vec(),
                                         sample_count);

                // Populate Base struct with known info tags
                if longread {
                    // get relevant flag for SVIM vcf on long read samples
                    base.variant = variant.clone();
                    base.depth[sample_idx] = match record.format(b"AD").integer() {
                        Ok(val) => {
                            if val[0][1] >= 0 {
                                val[0][1]
                            } else {
                                match record.info(b"SUPPORT").integer() {
                                    Ok(val) => {
                                        match val {
                                            Some(dep) => dep[0],
                                            _ => 0,
                                        }
                                    },
                                    _ => 0,
                                }
                            }
                        },
                        _ => {
                            0
                        }
                    };

//                    base.totaldepth[sample_idx] = match record.format(b"DP").integer() {
//                        Ok(val) => {
//                            if val[0][0] >= 0 {
//                                val[0][0]
//                            } else {
//                                base.depth[sample_idx]
//                            }
//                        },
//                        _ => base.depth[sample_idx],
//                    };
//                    let refr_depth = std::cmp::max(0, base.totaldepth[sample_idx] - base.depth[sample_idx]);
//                    base.af[sample_idx] = base.depth[sample_idx] as f64 / base.totaldepth[sample_idx] as f64;
//                    base.freq[sample_idx] = base.af[sample_idx];
                    let reads = record.info(b"READS").string().unwrap().unwrap().iter().map(|read| read.to_vec()).collect::<HashSet<Vec<u8>>>();
                    base.reads.par_extend(reads);
                    if refr_base_empty {
//                        let mut refr_base = Base::new(record.rid().unwrap(),
//                                                      record.pos(), record.alleles()[0].to_vec(), sample_count);
//                        refr_base.af[sample_idx] = base.af[sample_idx];
//                        refr_base.totaldepth[sample_idx] = match record.format(b"DP").integer() {
//                            Ok(val) => {
//                                if val[0][0] >= 0 {
//                                    val[0][0]
//                                } else {
//                                    base.depth[sample_idx]
//                                }
//                            },
//                            _ => base.depth[sample_idx],
//                        };
//                        refr_base.freq[sample_idx] = 1. - base.af[sample_idx];
//                        refr_base.depth[sample_idx] = refr_depth;
//                        bases.push(refr_base);
//                        refr_base_empty = false;
                    }
                } else {
                    // Get relevant flag from freebayes output on short read samples
                    base.variant = variant.clone();
                    base.filters[sample_idx] = filter_hash.clone();
//                    base.totaldepth[sample_idx] = record.info(b"DP").integer().unwrap().unwrap()[0];
                    base.baseq[sample_idx] = record.info(b"QA").integer().unwrap().unwrap()[0];
                    base.depth[sample_idx] = record.info(b"AO").integer().unwrap().unwrap()[0] as i32;
//                    base.referencedepth[sample_idx] = record.info(b"RO").integer().unwrap().unwrap()[0] as i32;

                    base.freq[sample_idx] = base.depth[sample_idx] as f64 / base.totaldepth[sample_idx] as f64;

                    if refr_base_empty {
//                        let mut refr_base = Base::new(record.rid().unwrap(),
//                                                      record.pos(), record.alleles()[0].to_vec(), sample_count);
//                        refr_base.totaldepth[sample_idx] = record.info(b"DP").integer().unwrap().unwrap()[0];
//                        refr_base.baseq[sample_idx] = record.info(b"QR").integer().unwrap().unwrap()[0];
//                        refr_base.depth[sample_idx] = record.info(b"RO").integer().unwrap().unwrap()[0] as i32;
//                        refr_base.freq[sample_idx] = refr_base.depth[sample_idx] as f64 / refr_base.totaldepth[sample_idx] as f64;
//
//                        bases.push(refr_base);
//                        refr_base_empty = false;
                    }
                };

                bases.push(base);
            };
            Some(bases)
        } else {
            None
        }
    }

    pub fn assign_read(&mut self, read_id: Vec<u8>) {
        self.reads.insert(read_id);
    }
}

/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    omit_snvs: bool,
    omit_indels: bool,
    indel_len_range: Option<Range<u32>>,
) -> Vec<Variant> {
    let pos = record.pos();
    let svlens = match record.info(b"SVLEN").integer() {
        Ok(Some(svlens)) => Some(
            svlens
                .into_iter()
                .map(|l| {
                    if !l.is_missing() {
                        Some(l.abs() as u32)
                    } else {
                        None
                    }
                })
                .collect_vec(),
        ),
        _ => None,
    };
    let end = match record.info(b"END").integer() {
        Ok(Some(end)) => {
            let end = end[0] as u32 - 1;
            Some(end)
        }
        _ => None,
    };
    // TODO avoid cloning svtype
    let svtype = match record.info(b"SVTYPE").string() {
        Ok(Some(svtype)) => Some(svtype[0].to_owned()),
        _ => None,
    };

    // check if len is within the given range
    let is_valid_len = |svlen| {
        if let Some(ref len_range) = indel_len_range {
            // TODO replace with Range::contains once stabilized
            if svlen < len_range.start || svlen >= len_range.end {
                return false;
            }
        }
        true
    };

    let is_valid_insertion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<INS>"
            || (ref_allele.len() < alt_allele.len()
            && ref_allele == &alt_allele[..ref_allele.len()])
    };

    let is_valid_deletion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<DEL>"
            || (ref_allele.len() > alt_allele.len()
            && &ref_allele[..alt_allele.len()] == alt_allele)
    };

    let variants = if let Some(svtype) = svtype {
        vec![if omit_indels {
            Variant::None
        } else if svtype == b"INS" {
            // get sequence
            let alleles = record.alleles();
            if alleles.len() > 2 {
                panic!("SVTYPE=INS but more than one ALT allele".to_owned());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];

            if alt_allele == b"<INS>" {
                // don't support insertions without exact sequence
                Variant::None
            } else {
                let len = alt_allele.len() - ref_allele.len();

                if is_valid_insertion_alleles(ref_allele, alt_allele) && is_valid_len(len as u32) {
                    Variant::Insertion(
                        alt_allele[ref_allele.len()..].to_owned(),
                    )
                } else {
                    Variant::None
                }
            }
        } else if svtype == b"DEL" {
            let svlen = match (svlens, end) {
                (Some(ref svlens), _) if svlens[0].is_some() => svlens[0].unwrap(),
                (None, Some(end)) => end - (pos as u32 + 1), // pos is pointing to the allele before the DEL
                _ => {
                    panic!("SVLEN or END".to_owned());
                }
            };
            if svlen == 0 {
                panic!("Absolute value of SVLEN or END - POS must be greater than zero."
                    .to_owned());
            }
            let alleles = record.alleles();
            if alleles.len() > 2 {
                panic!("SVTYPE=DEL but more than one ALT allele".to_owned());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];

            if alt_allele == b"<DEL>" || is_valid_deletion_alleles(ref_allele, alt_allele) {
                if is_valid_len(svlen) {
                    Variant::Deletion(svlen)
                } else {
                    Variant::None
                }
            } else {
                Variant::None
            }
        } else {
            Variant::None
        }]
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];
        let mut variant_vec = vec!();
        alleles
            .iter()
            .skip(1)
            .enumerate()
            .for_each(|(i, alt_allele)| {
                if alt_allele == b"<*>" {
                    // dummy non-ref allele, signifying potential homozygous reference site
                    if omit_snvs {
                        variant_vec.push(Variant::None)
                    } else {
                        variant_vec.push(Variant::None)
                    }
                } else if alt_allele == b"<DEL>" {
                    if let Some(ref svlens) = svlens {
                        if let Some(svlen) = svlens[i] {
                            variant_vec.push(Variant::Deletion(svlen))
                        } else {
                            // TODO fail with an error in this case
                            variant_vec.push(Variant::None)
                        }
                    } else {
                        // TODO fail with an error in this case
                        variant_vec.push(Variant::None)
                    }
                } else if alt_allele[0] == b'<' {
                    // TODO Catch <DUP> structural variants here
                    // skip any other special alleles
                    variant_vec.push(Variant::None)
                } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                    // SNV
                    if omit_snvs {
                        variant_vec.push(Variant::None)
                    } else if alt_allele == &ref_allele {
                        variant_vec.push(Variant::None)
                    } else {
                        variant_vec.push(Variant::SNV(alt_allele[0]))
                    }
                } else if alt_allele.len() == ref_allele.len() {
                    // MNV
                    variant_vec.push(Variant::MNV(alt_allele.to_vec()))
                } else {
                    let indel_len =
                        (alt_allele.len() as i32 - ref_allele.len() as i32).abs() as u32;
                    // TODO fix position if variant is like this: cttt -> ct

                    if omit_indels || !is_valid_len(indel_len) {
                        variant_vec.push(Variant::MNV(alt_allele.to_vec()))
                    } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                        variant_vec.push(Variant::Deletion(
                            (ref_allele.len() - alt_allele.len()) as u32,
                        ))
                    } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                        variant_vec.push(Variant::Insertion(
                            alt_allele[ref_allele.len()..].to_owned(),
                        ))
                    } else {
                        variant_vec.push(Variant::MNV(alt_allele.to_vec()))
                    }
                }
            });
        variant_vec
    };

    variants
}

#[cfg(test)]
mod tests {
    use super::*;


}
