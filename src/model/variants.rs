use bio::stats::LogProb;
use itertools::Itertools;
use ordered_float::{NotNan, OrderedFloat};
use rust_htslib::{bcf, bcf::record::Numeric};
use std::collections::HashSet;
use std::fmt::Debug;
use std::ops::Range;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};
use utils::ReadType;

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
    #[strum(serialize = "INV")]
    Inversion(Option<Range<u32>>),
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
            "INV" => VariantType::Inversion(None),
            "SNV" => VariantType::SNV,
            "REF" => VariantType::None,
            _ => panic!("bug: given string does not describe a valid variant type"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub enum SVType {
    DUP,
    INV,
}

#[derive(Clone, Copy, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub struct SV {
    sv: SVType,
    len: u32,
    start: u32,
    end: u32,
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq, Sized)]
pub enum Variant {
    Deletion(u32),
    Insertion(Vec<u8>),
    Inversion(Vec<u8>),
    SNV(u8),
    MNV(Vec<u8>),
    SV(SV),
    None,
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq, Sized)]
pub struct Allele {
    variant: Variant,
    reference: bool
}

impl Allele {

    const SINGLE_BREAKEND_INDICATOR: char = '.';
    const BREAKEND_EXTENDING_RIGHT: char = '[';
    const BREAKEND_EXTENDING_LEFT: char = ']';
    const SYMBOLIC_ALLELE_START: char = '<';
    const SYMBOLIC_ALLELE_END: char = '>';

    const NO_CALL: char = '.';
    const SPAND_DEL: char = '*';

    const NON_REF_STRING: String = "<NON_REF>".to_string();
    pub const NON_REF_ALLELE: Allele = Allele::new(
        Variant::MNV(Allele::NON_REF_STRING.into_bytes()),
        false
    );

    pub fn new(variant: Variant, reference: bool) -> Allele {
        Allele {
            variant,
            reference,
        }
    }

    pub fn create_fake_alleles() -> Vec<Allele> {
        let mut alleles = vec![Allele::fake(true), Allele::fake(false)];

        return alleles
    }

    pub fn fake(reference: bool) -> Allele {
        Allele {
            variant: Variant::None,
            reference,
        }
    }

    pub fn no_call() -> Allele {
        Allele {
            variant: Variant::SNV(Allele::NO_CALL as u8),
            reference: false,
        }
    }

    pub fn is_reference(&self) -> bool {
        self.reference
    }

    pub fn unwrap(possible_allele: Option<&Allele>) -> Allele {
        let a = match possible_allele {
            Some(a) => {
                *a.clone()
            },
            _ => Allele::fake(false)
        };
        a
    }

    pub fn is_no_call(&self) -> bool {
        match &self.variant {
            Variant::SNV(snp) => {
                snp == &(Allele::NO_CALL as u8)
            },
            _ => false
        }
    }

    pub fn is_called(&self) -> bool {
        match &self.variant {
            Variant::SNV(snp) => {
                snp != &(Allele::NO_CALL as u8)
            },
            _ => true
        }
    }

    pub fn is_symbolic(&self) -> bool {
        self.variant.is_symbolic()
    }

    pub fn length(&self) -> usize {
        self.variant.len()
    }

    pub fn variant(&self) -> &Variant {
        &self.variant
    }

    pub fn contains(&self, variant: &Variant) -> bool {
        self.variant == variant
    }

    pub fn is_del(&self) -> bool {
        match self.variant {
            &Variant::Deletion(_) => true,
            _ => false,
        }
    }
}


#[allow(unused)]
impl Variant {

    const SINGLE_BREAKEND_INDICATOR: u8 = '.' as u8;
    const BREAKEND_EXTENDING_RIGHT: u8 = '[' as u8;
    const BREAKEND_EXTENDING_LEFT: u8 = ']' as u8;
    const SYMBOLIC_ALLELE_START: u8 = '<' as u8;
    const SYMBOLIC_ALLELE_END: u8 = '>' as u8;

    const NO_CALL: u8 = '.' as u8;
    const SPAND_DEL: u8 = '*' as u8;

    pub fn is_symbolic(&self) -> bool {
        match self {
            &Variant::Inversion(var)
            | &Variant::Insertion(var)
            | &Variant::MNV(var) => {
                if var.contains(Variant::SINGLE_BREAKEND_INDICATOR)
                    | var.contains(Variant::BREAKEND_EXTENDING_LEFT)
                    | var.contains(Variant::BREAKEND_EXTENDING_RIGHT)
                    | var.contains(Variant::SYMBOLIC_ALLELE_START)
                    | var.contains(Variant::SYMBOLIC_ALLELE_END)
                    | var.contains(Variant::NO_CALL) {
                    true
                } else {
                    false
                }
            },
            &Variant::SNV(var) => {
                if (var == Variant::SINGLE_BREAKEND_INDICATOR)
                    | (var == Variant::BREAKEND_EXTENDING_LEFT)
                    | (var == Variant::BREAKEND_EXTENDING_RIGHT)
                    | (var == Variant::SYMBOLIC_ALLELE_START)
                    | (var == Variant::SYMBOLIC_ALLELE_END)
                    | (var == Variant::NO_CALL) {
                    true
                } else {
                    false
                }
            },
            &Variant::Deletion(_)
            | &Variant::SV(_) => {
                false
            }
            &Variant::None => true,

        }
    }

    pub fn has_fragment_evidence(&self) -> bool {
        match self {
            &Variant::Deletion(_) => true,
            &Variant::Insertion(_) => true,
            &Variant::Inversion(_) => true,
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
            &Variant::Inversion(_) => false,
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
            (&Variant::Inversion(_), &VariantType::Inversion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            }
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::Inversion(_), &VariantType::Inversion(None)) => true,
            (&Variant::SNV(_), &VariantType::SNV) => true,
            (&Variant::MNV(_), &VariantType::MNV) => true,
            (&Variant::None, &VariantType::None) => true,
            _ => false,
        }
    }

    pub fn end(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length,
            &Variant::Insertion(_) | &Variant::Inversion(_) => start + 1, // end of insertion is the next regular base
            &Variant::SV(sv) => sv.end,
            &Variant::SNV(_) | &Variant::None => start,
            &Variant::MNV(ref alt) => start + alt.len() as u32,
        }
    }

    pub fn centerpoint(&self, start: u32) -> u32 {
        match self {
            &Variant::Deletion(length) => start + length / 2,
            &Variant::Insertion(_) | &Variant::Inversion(_) => start, // end of insertion is the next regular base
            &Variant::SV(sv) => (sv.start + sv.len) / 2,
            &Variant::SNV(_) | &Variant::None => start,
            &Variant::MNV(ref alt) => start + alt.len() as u32 / 2,
        }
    }

    pub fn len(&self) -> u32 {
        match self {
            &Variant::Deletion(l) => l,
            &Variant::Insertion(ref s) | &Variant::Inversion(ref s) => s.len() as u32,
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
    // The QUAL values across samples
    pub quals: Vec<f64>,
    // Depth of this variant as decided by variant caller. Only includes good quality reads
    // This values tends to be inconsistent so opt to use lorikeets depth value
    pub depth: Vec<i32>,
    // Depth for this variant as decided by lorikeet. Includes low quality reads
    pub truedepth: Vec<i32>,
    // Depth at this location as decided by CoverM
    pub totaldepth: Vec<i32>,
    // Depth of the reference allele as decided by lorikeet. Includes low quality reds
    pub referencedepth: Vec<i32>,
    // Frequency of variant
    pub freq: Vec<f64>,
    // Read ids assigned to variant
    pub reads: HashSet<Vec<u8>>,
    // CLR transformed relative abundances
    pub rel_abunds: Vec<f64>,
    // Genotypes assigned to variant
    pub genotypes: HashSet<i32>,
}

#[allow(unused)]
impl Base {
    /// Update depth as calculated by lorikeet and the reference depth
    pub fn add_depth(&mut self, sample_idx: usize, d: i32, refr_depth: i32) {
        if self.totaldepth[sample_idx] == 0 {
            self.totaldepth[sample_idx] = d;
            match self.variant {
                _ => {
                    self.referencedepth[sample_idx] = refr_depth;
                }
            }
        }
    }

    pub fn combine_sample(&mut self, other: &Base, sample_idx: usize, total_depth: i32) {
        if &self != &other {
            self.quals[sample_idx] = other.quals[sample_idx];
            self.depth[sample_idx] = other.depth[sample_idx];
            self.truedepth[sample_idx] = other.truedepth[sample_idx];
            self.totaldepth[sample_idx] = total_depth;
            self.freq[sample_idx] = other.freq[sample_idx];
        } else {
            self.totaldepth[sample_idx] = total_depth;
        }
    }

    pub fn new(tid: u32, pos: i64, sample_count: usize, refr: Vec<u8>) -> Base {
        Base {
            tid,
            pos,
            refr,
            variant: Variant::None,
            quals: vec![0.; sample_count],
            depth: vec![0; sample_count],
            truedepth: vec![0; sample_count],
            totaldepth: vec![0; sample_count],
            referencedepth: vec![0; sample_count],
            freq: vec![0.; sample_count],
            rel_abunds: vec![0.; sample_count],
            reads: HashSet::new(),
            genotypes: HashSet::new(),
        }
    }

    pub fn from_vcf_record(
        record: &mut bcf::Record,
        sample_count: usize,
        sample_idx: usize,
        readtype: &ReadType,
        min_qual: f32,
    ) -> Option<Vec<Base>> {
        if record.qual() > min_qual {
            let variants = collect_variants(record, false, false, None);
            if variants.len() > 0 {
                // Separate filters into hashset of filter struct
                let mut filter_hash = HashSet::new();
                {
                    let filters = record.filters();
                    let header = record.header();
                    for filter in filters {
                        filter_hash.insert(Filter::from_result(std::str::from_utf8(
                            &header.id_to_name(filter)[..],
                        )));
                    }
                }
                let mut bases = vec![];
                let mut refr_base_empty = true;
                for (idx, variant) in variants.iter().enumerate() {
                    // Get elements from record
                    let mut base = Base::new(
                        record.rid().unwrap(),
                        record.pos(),
                        sample_count,
                        record.alleles()[0].to_vec(),
                    );
                    base.quals[sample_idx] = record.qual() as f64;

                    // Populate Base struct with known info tags
                    match readtype {
                        &ReadType::Long | &ReadType::Assembly => {
                            // get relevant flag for SVIM vcf on long read samples
                            base.variant = variant.clone();
                            base.depth[sample_idx] = match record.format(b"AD").integer() {
                                Ok(val) => {
                                    if val[0][1] >= 0 {
                                        val[0][1]
                                    } else {
                                        match record.info(b"SUPPORT").integer() {
                                            Ok(val) => match val {
                                                Some(dep) => dep[0],
                                                _ => 0,
                                            },
                                            _ => 0,
                                        }
                                    }
                                }
                                _ => 0,
                            };

                            base.truedepth[sample_idx] = match record.format(b"AD").integer() {
                                Ok(val) => {
                                    if val[0][1] >= 0 {
                                        val[0][1]
                                    } else {
                                        match record.info(b"SUPPORT").integer() {
                                            Ok(val) => match val {
                                                Some(dep) => dep[0],
                                                _ => 0,
                                            },
                                            _ => 0,
                                        }
                                    }
                                }
                                _ => 0,
                            };

                            let refr_depth = std::cmp::max(
                                0,
                                base.totaldepth[sample_idx] - base.depth[sample_idx],
                            );
                            //                    base.af[sample_idx] = base.depth[sample_idx] as f64 / base.totaldepth[sample_idx] as f64;
                            //                    base.freq[sample_idx] = base.af[sample_idx];
                            let reads = record
                                .info(b"READS")
                                .string()
                                .unwrap()
                                .unwrap()
                                .iter()
                                .map(|read| read.to_vec())
                                .collect::<HashSet<Vec<u8>>>();
                            base.reads.par_extend(reads);
                            if refr_base_empty {
                                let mut refr_base = Base::new(
                                    record.rid().unwrap(),
                                    record.pos(),
                                    sample_count,
                                    record.alleles()[0].to_vec(),
                                );

                                refr_base.depth[sample_idx] = refr_depth;
                                bases.push(refr_base);
                                refr_base_empty = false;
                            }
                        }
                        &ReadType::Short => {
                            // Get relevant flag from freebayes output on short read samples
                            // let allele_depths = record.format(b"AD").integer().unwrap().clone();
                            base.variant = variant.clone();
                            //                    base.totaldepth[sample_idx] = record.info(b"DP").integer().unwrap().unwrap()[0];
                            // base.baseq[sample_idx] = record.info(b"QA").integer().unwrap().unwrap()[0];
                            base.depth[sample_idx] = match record.format(b"AD").integer() {
                                Ok(val) => val[0][1],
                                _ => {
                                    debug!("No AD format");

                                    0
                                }
                            };
                            //                    base.referencedepth[sample_idx] = record.info(b"RO").integer().unwrap().unwrap()[0] as i32;

                            //                        base.freq[sample_idx] = base.depth[sample_idx] as f64 / base.totaldepth[sample_idx] as f64;

                            if refr_base_empty {
                                let mut refr_base = Base::new(
                                    record.rid().unwrap(),
                                    record.pos(),
                                    sample_count,
                                    record.alleles()[0].to_vec(),
                                );
                                //                        refr_base.totaldepth[sample_idx] = record.info(b"DP").integer().unwrap().unwrap()[0];
                                // refr_base.baseq[sample_idx] =
                                //     record.info(b"QR").integer().unwrap().unwrap()[0];
                                refr_base.depth[sample_idx] = match record.format(b"AD").integer() {
                                    Ok(val) => val[0][0],
                                    _ => 0,
                                };
                                //                        refr_base.freq[sample_idx] = refr_base.depth[sample_idx] as f64 / refr_base.totaldepth[sample_idx] as f64;
                                //
                                bases.push(refr_base);
                                refr_base_empty = false;
                            }
                        }
                    }
                    debug!("Variant {:?}", base.variant);
                    bases.push(base);
                }
                Some(bases)
            } else {
                None
            }
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
        // Gets value from SVLEN tag in VCF record
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

    let is_valid_inversion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<INV>" || (ref_allele.len() == alt_allele.len())
    };

    let is_valid_mnv = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<MNV>" || (ref_allele.len() == alt_allele.len())
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
                    Variant::Insertion(alt_allele[ref_allele.len()..].to_owned())
                } else {
                    Variant::None
                }
            }
        } else if svtype == b"INV" {
            // get sequence
            let alleles = record.alleles();
            if alleles.len() > 2 {
                panic!("SVTYPE=INS but more than one ALT allele".to_owned());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];
            debug!("Inversion being collected...");
            if alt_allele == b"<INV>" {
                // don't support inversions without exact sequence
                Variant::None
            } else {
                let len = alt_allele.len();

                if is_valid_inversion_alleles(ref_allele, alt_allele) && is_valid_len(len as u32) {
                    Variant::Inversion(alt_allele.to_owned())
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
                panic!(
                    "Absolute value of SVLEN or END - POS must be greater than zero.".to_owned()
                );
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
        let mut variant_vec = vec![];
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
                } else {
                    let indel_len =
                        (alt_allele.len() as i32 - ref_allele.len() as i32).abs() as u32;
                    // TODO fix position if variant is like this: cttt -> ct

                    if (omit_indels || !is_valid_len(indel_len))
                        && is_valid_mnv(ref_allele, alt_allele)
                    {
                        // println!("MNV 1 {:?} {:?}", ref_allele, alt_allele);
                        variant_vec.push(Variant::MNV(alt_allele.to_vec()))
                    } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                        variant_vec.push(Variant::Deletion(
                            (ref_allele.len() - alt_allele.len()) as u32,
                        ))
                    } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                        variant_vec.push(Variant::Insertion(
                            alt_allele[ref_allele.len()..].to_owned(),
                        ))
                    } else if is_valid_mnv(ref_allele, alt_allele) {
                        // println!("MNV 2 {:?} {:?}", ref_allele, alt_allele);
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
