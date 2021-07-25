use haplotype::event_map::EventMap;
use model::byte_array_allele::ByteArrayAllele;
use ordered_float::{NotNan, OrderedFloat};
use rust_htslib::bam::record::{Cigar, CigarString};
use std::hash::{Hash, Hasher};
use utils::simple_interval::Locatable;

// lazy_static! {
//     pub static ref SIZE_AND_BASE_ORDER: Then<Extract<Fn(&Haplotype<Locatable>)>>
// }
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Haplotype<L: Locatable> {
    pub(crate) allele: ByteArrayAllele,
    pub(crate) genome_location: Option<L>,
    pub(crate) event_map: Option<EventMap<L>>,
    pub(crate) cigar: CigarString,
    pub(crate) alignment_start_hap_wrt_ref: usize,
    pub(crate) score: OrderedFloat<f64>,
    // debug information for tracking kmer sizes used in graph construction for debug output
    pub(crate) kmer_size: usize,
}

impl<L: Locatable> Haplotype<L> {
    /**
     * Main constructor
     *
     * @param bases a non-null array of bases
     * @param isRef is this the reference haplotype?
     */
    pub fn new(bases: &[u8], is_ref: bool) -> Haplotype<L> {
        Haplotype {
            allele: ByteArrayAllele::new(bases, is_ref),
            genome_location: None,
            event_map: None,
            cigar: CigarString::from(vec![Cigar::Match(0)]),
            alignment_start_hap_wrt_ref: 0,
            score: OrderedFloat(std::f64::MIN),
            kmer_size: 0,
        }
    }

    pub fn set_alignment_start_hap_wrt_ref(&mut self, value: usize) {
        self.alignment_start_hap_wrt_ref = value
    }

    pub fn get_bases(&self) -> &Vec<u8> {
        &self.allele.bases
    }

    pub fn set_cigar(&mut self, cigar_string: Vec<Cigar>) {
        self.cigar = CigarString::from(cigar_string)
    }

    pub fn set_genome_location(&mut self, genome_location: L) {
        self.genome_location = Some(genome_location)
    }

    pub fn len(&self) -> usize {
        self.allele.len()
    }
}

impl<L: Locatable> Hash for Haplotype<L> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cigar.hash(state);
        self.allele.hash(state);
        self.genome_location.hash(state);
        self.score.hash(state);
        self.kmer_size.hash(state);
    }
}
