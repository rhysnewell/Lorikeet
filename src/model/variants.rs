use bio::stats::LogProb;
use itertools::Itertools;
use ordered_float::NotNan;
use std::collections::HashSet;
use std::fmt::Debug;
use std::ops::Range;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};

use rayon::prelude::*;
use model::variant_context::VariantContext;
use utils::vcf_constants::VCFConstants;
use rayon::prelude::*;

pub type AlleleFreq = NotNan<f64>;

lazy_static! {
    static ref NON_REF_STRING: String = "<NON_REF>".to_string();
    pub static ref NON_REF_ALLELE: Allele = Allele::new(
        Variant::MNV((&(*NON_REF_STRING)).clone().into_bytes().to_vec()),
        false
    );
}

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

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub enum Variant {
    Deletion(u32),
    Insertion(Vec<u8>),
    Inversion(Vec<u8>),
    SNV(u8),
    MNV(Vec<u8>),
    SV(SV),
    None,
}

#[derive(Clone, Debug, PartialEq, Ord, PartialOrd, Hash, Eq)]
pub struct Allele {
    variant: Variant,
    reference: bool
}


/**
 * Immutable representation of an allele.
 *<p>
 * Types of alleles:
 *</p>
 *<pre>
 Ref: a t C g a // C is the reference base
 : a t G g a // C base is a G in some individuals
 : a t - g a // C base is deleted w.r.t. the reference
 : a t CAg a // A base is inserted w.r.t. the reference sequence
 </pre>
 *<p> In these cases, where are the alleles?</p>
 *<ul>
 * <li>SNP polymorphism of C/G  -&gt; { C , G } -&gt; C is the reference allele</li>
 * <li>1 base deletion of C     -&gt; { tC , t } -&gt; C is the reference allele and we include the preceding reference base (null alleles are not allowed)</li>
 * <li>1 base insertion of A    -&gt; { C ; CA } -&gt; C is the reference allele (because null alleles are not allowed)</li>
 *</ul>
 *<p>
 * Suppose I see a the following in the population:
 *</p>
 *<pre>
 Ref: a t C g a // C is the reference base
 : a t G g a // C base is a G in some individuals
 : a t - g a // C base is deleted w.r.t. the reference
 </pre>
 * <p>
 * How do I represent this?  There are three segregating alleles:
 * </p>
 *<blockquote>
 *  { C , G , - }
 *</blockquote>
 *<p>and these are represented as:</p>
 *<blockquote>
 *  { tC, tG, t }
 *</blockquote>
 *<p>
 * Now suppose I have this more complex example:
 </p>
 <pre>
 Ref: a t C g a // C is the reference base
 : a t - g a
 : a t - - a
 : a t CAg a
 </pre>
 * <p>
 * There are actually four segregating alleles:
 * </p>
 *<blockquote>
 *   { Cg , -g, --, and CAg } over bases 2-4
 *</blockquote>
 *<p>   represented as:</p>
 *<blockquote>
 *   { tCg, tg, t, tCAg }
 *</blockquote>
 *<p>
 * Critically, it should be possible to apply an allele to a reference sequence to create the
 * correct haplotype sequence:</p>
 *<blockquote>
 * Allele + reference =&gt; haplotype
 *</blockquote>
 *<p>
 * For convenience, we are going to create Alleles where the GenomeLoc of the allele is stored outside of the
 * Allele object itself.  So there's an idea of an A/C polymorphism independent of it's surrounding context.
 *
 * Given list of alleles it's possible to determine the "type" of the variation
 </p>
 <pre>
 A / C @ loc =&gt; SNP
 - / A =&gt; INDEL
 </pre>
 * <p>
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion.
 * </p>
 * <p>
 * Alelle also supports is concept of a NO_CALL allele.  This Allele represents a haplotype that couldn't be
 * determined. This is usually represented by a '.' allele.
 * </p>
 * <p>
 * Note that Alleles store all bases as bytes, in **UPPER CASE**.  So 'atc' == 'ATC' from the perspective of an
 * Allele.
 * </p>
 * @author gatk_team.
 */
impl Allele {

    const SINGLE_BREAKEND_INDICATOR: char = '.';
    const BREAKEND_EXTENDING_RIGHT: char = '[';
    const BREAKEND_EXTENDING_LEFT: char = ']';
    const SYMBOLIC_ALLELE_START: char = '<';
    const SYMBOLIC_ALLELE_END: char = '>';

    const NO_CALL: char = '.';
    const SPAND_DEL: char = '*';


    /**
     * Create a new Allele that includes bases and if tagged as the reference allele if isRef == true.  If bases
     * == '-', a Null allele is created.  If bases ==  '.', a no call Allele is created. If bases ==  '*', a spanning deletions Allele is created.
     *
     * @param bases the DNA sequence of this variation, '-', '.', or '*'
     * @param isRef should we make this a reference allele?
     * @throws IllegalArgumentException if bases contains illegal characters or is otherwise malformated
     */
    pub fn create(bases: &[u8], is_ref: bool) -> Allele {
        if bases.len() == 1 {
            let base = bases[0].to_ascii_uppercase();
            return Allele::new(Variant::SNV(base), is_ref)
        } else {

        }
    }

    pub fn new(variant: Variant, reference: bool) -> Allele {
        Allele {
            variant,
            reference,
        }
    }

    pub fn create_fake_alleles() -> Vec<Allele> {
        let alleles = vec![Allele::fake(true), Allele::fake(false)];

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
                a.clone()
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
        self.variant.len() as usize
    }

    pub fn variant(&self) -> &Variant {
        &self.variant
    }

    pub fn contains(&self, variant: &Variant) -> bool {
        &self.variant == variant
    }

    pub fn is_del(&self) -> bool {
        match self.variant {
            Variant::Deletion(_) => true,
            _ => false,
        }
    }

    pub fn would_be_null_allele(bases: &[u8]) -> bool {
        return bases.len() == 1 && bases[0] as char == VCFConstants::NULL_ALLELE || bases.len() == 0
    }

    pub fn would_be_no_call_allele(bases: &[u8]) -> bool {
        return bases.len() == 1 && bases[0] as char == VCFConstants::NO_CALL_ALLELE
    }

    pub fn would_be_star_allele(bases: &[u8]) -> bool {
        return bases.len() == 1 && bases[0] as char == VCFConstants::SPANNING_DELETION_ALLELE
    }

    pub fn would_be_symbolic_allele(bases: &[u8]) -> bool {
        if bases.len() <= 1 {
            return false
        } else {
            return bases[0] as char == Self::SYMBOLIC_ALLELE_START || bases[bases.len() - 1] == Self::SYMBOLIC_ALLELE_END ||
                Self::would_be_breakpoint(bases) || Self::would_be_single_breakend(bases)
        }
    }

    pub fn would_be_breakpoint(bases: &[u8]) -> bool {
        if bases.len() <= 1 {
            return false
        }
        return bases.iter().par_bridge().any(|base| base as char == Self::BREAKEND_EXTENDING_LEFT || base as char == Self::BREAKEND_EXTENDING_RIGHT);
    }

    pub fn would_be_single_breakend(bases: &[u8]) -> bool {
        if bases.len() <= 1 {
            return false
        } else {
            return bases[0] == Self::SINGLE_BREAKEND_INDICATOR || bases[bases.len() - 1] == Self::SINGLE_BREAKEND_INDICATOR
        }
    }

    pub fn acceptable_allele_bases(bases: &[u8], is_ref: bool) -> bool {
        if Self::would_be_null_allele(bases) {
            return false
        } else if Self::would_be_no_call_allele(bases) || Self::would_be_symbolic_allele(bases) {
            return true
        } else if Self::would_be_star_allele(bases) {
            return !is_ref
        } else {
            // return true if there are any unacceptable bases, so take conjugate value
            !bases.iter().par_bridge().any(|base| {
                let base = base as char;
                match base {
                    'A' | 'C' | 'T' | 'G' |
                    'a' | 'c' | 't' | 'g' |
                    'N' | 'n' => false,
                    _ => true
                }
            })
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
            Variant::Inversion(var)
            | Variant::Insertion(var)
            | Variant::MNV(var) => {
                if var.contains(&Variant::SINGLE_BREAKEND_INDICATOR)
                    | var.contains(&Variant::BREAKEND_EXTENDING_LEFT)
                    | var.contains(&Variant::BREAKEND_EXTENDING_RIGHT)
                    | var.contains(&Variant::SYMBOLIC_ALLELE_START)
                    | var.contains(&Variant::SYMBOLIC_ALLELE_END)
                    | var.contains(&Variant::NO_CALL) {
                    true
                } else {
                    false
                }
            },
            Variant::SNV(var) => {
                if (var == &Variant::SINGLE_BREAKEND_INDICATOR)
                    | (var == &Variant::BREAKEND_EXTENDING_LEFT)
                    | (var == &Variant::BREAKEND_EXTENDING_RIGHT)
                    | (var == &Variant::SYMBOLIC_ALLELE_START)
                    | (var == &Variant::SYMBOLIC_ALLELE_END)
                    | (var == &Variant::NO_CALL) {
                    true
                } else {
                    false
                }
            },
            Variant::Deletion(_)
            | Variant::SV(_) => {
                false
            },
            Variant::None => true,

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

#[cfg(test)]
mod tests {
    use super::*;
}
