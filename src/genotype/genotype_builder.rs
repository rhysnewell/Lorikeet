use genotype::genotype_likelihoods::GenotypeLikelihoods;
use model::byte_array_allele::ByteArrayAllele;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

lazy_static! {
    static ref HAPLOID_NO_CALL: Vec<ByteArrayAllele> = vec![ByteArrayAllele::fake(false)];
    static ref DIPLOID_NO_CALL: Vec<ByteArrayAllele> = vec![ByteArrayAllele::fake(false); 2];
}

#[derive(Debug, Clone)]
pub enum GenotypeAssignmentMethod {
    BestMatchToOriginal,
    DoNotAssignGenotypes,
    SetToNoCall,
    SetToNoCallNoAnnotations,
    UsePLsToAssign,
    UsePosteriorProbabilities,
}

impl GenotypeAssignmentMethod {
    pub fn from_args(args: &clap::ArgMatches) -> GenotypeAssignmentMethod {
        match args.value_of("genotype-assignment-method").unwrap() {
            "UsePLsToAssign" => GenotypeAssignmentMethod::UsePLsToAssign,
            "UsePosteriorProbabilities" => GenotypeAssignmentMethod::UsePosteriorProbabilities,
            "BestMatchToOriginal" => GenotypeAssignmentMethod::BestMatchToOriginal,
            _ => GenotypeAssignmentMethod::DoNotAssignGenotypes,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum AttributeObject {
    f64(f64),
    Vecf64(Vec<f64>),
    String(String),
    UnsizedInteger(usize),
    VecUnsize(Vec<usize>),
    VecU8(Vec<u8>),
    None,
}

impl AttributeObject {
    pub fn is_empty(&self) -> bool {
        match self {
            Self::Vecf64(vec) => vec.is_empty(),
            Self::String(string) => string.is_empty(),
            Self::None => true,
            _ => false,
        }
    }

    pub fn to_ref_vec(&self) -> &Vec<f64> {
        match self {
            Self::Vecf64(vec) => &vec,
            _ => panic!("Attrribute object can't be cast to vec"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Genotype {
    pub ploidy: usize,
    pub pl: Vec<i64>,
    pub alleles: Vec<ByteArrayAllele>,
    pub ad: Vec<i64>,
    pub dp: i64,
    pub gq: i64,
    pub is_phased: bool,
    pub sample_name: String,
    pub attributes: HashMap<String, AttributeObject>,
    pub genotype_type: Option<GenotypeType>,
}

impl Eq for Genotype {}

impl PartialEq for Genotype {
    fn eq(&self, other: &Self) -> bool {
        self.ploidy == other.ploidy
            && self.alleles == other.alleles
            && self.ad == other.ad
            && self.dp == other.dp
            && self.gq == other.gq
            && self.is_phased == other.is_phased
            && self.sample_name == other.sample_name
    }
}

impl Hash for Genotype {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ploidy.hash(state);
        self.alleles.hash(state);
        self.ad.hash(state);
        self.dp.hash(state);
        self.gq.hash(state);
        self.is_phased.hash(state);
        self.sample_name.hash(state);
    }
}

impl Genotype {
    pub fn build(default_ploidy: usize, likelihoods: Vec<f64>, sample_name: String) -> Genotype {
        Genotype {
            ploidy: default_ploidy,
            alleles: Vec::with_capacity(likelihoods.len()),
            ad: Vec::with_capacity(likelihoods.len()),
            pl: GenotypeLikelihoods::from_log10_likelihoods(likelihoods).as_pls(),
            dp: -1,
            gq: -1,
            is_phased: false,
            sample_name,
            attributes: HashMap::new(),
            genotype_type: None,
        }
    }

    pub fn build_from_likelihoods(
        default_ploidy: usize,
        likelihoods: GenotypeLikelihoods,
        sample_name: String,
    ) -> Genotype {
        Genotype {
            ploidy: default_ploidy,
            alleles: Vec::with_capacity(likelihoods.len()),
            ad: Vec::with_capacity(likelihoods.len()),
            pl: likelihoods.as_pls(),
            dp: -1,
            gq: -1,
            is_phased: false,
            sample_name,
            attributes: HashMap::new(),
            genotype_type: None,
        }
    }

    pub fn build_from_alleles(alleles: Vec<ByteArrayAllele>, sample_name: String) -> Genotype {
        Genotype {
            ploidy: alleles.len(),
            pl: Vec::with_capacity(alleles.len()),
            dp: -1,
            gq: -1,
            ad: Vec::with_capacity(alleles.len()),
            is_phased: false,
            attributes: HashMap::new(),
            sample_name,
            alleles,
            genotype_type: None,
        }
    }

    pub fn get_ploidy(&self) -> usize {
        self.ploidy
    }

    pub fn get_likelihoods(&self) -> GenotypeLikelihoods {
        GenotypeLikelihoods::from_pls(self.pl.clone())
    }

    // pub fn get_likelihoods_mut(&mut self) -> &mut GenotypeLikelihoods {
    //     &mut self.pl
    // }

    /**
     * Compute how many likelihood elements are associated with the given number of alleles
     * Equivalent to asking in how many ways N non-negative integers can add up to P is S(N,P)
     * where P = ploidy (number of chromosomes copies) and N = total # of alleles.
     * Each chromosome can be in one single state (0,...,N-1) and there are P of them.
     * Naive solution would be to store N*P likelihoods, but this is not necessary because we can't distinguish chromosome states, but rather
     * only total number of alt allele counts in all chromosomes.
     *
     * For example, S(3,2) = 6: For alleles A,B,C, on a diploid organism we have six possible genotypes:
     * AA,AB,BB,AC,BC,CC.
     * Another way of expressing is with vector (#of A alleles, # of B alleles, # of C alleles)
     * which is then, for ordering above, (2,0,0), (1,1,0), (0,2,0), (1,1,0), (0,1,1), (0,0,2)
     * In general, for P=2 (regular biallelic), then S(N,2) = N*(N+1)/2
     *
     *
     * the result is calculated via a call to calcNumLikelihoods,
     * which uses the Apache Commons CombinatoricsUtils class
     * using the formula (numAlleles + ploidy - 1) choose ploidy
     *
     *   @param  numAlleles      Number of alleles (including ref)
     *   @param  ploidy          Ploidy, or number of chromosomes in set
     *   @return    Number of likelihood elements we need to hold.
     */
    pub fn num_likelihoods(&mut self, num_alleles: i64, ploidy: i64) -> usize {
        let result = self.get_likelihoods().num_likelihoods(num_alleles, ploidy);
        if result < 0 {
            0
        } else {
            result as usize
        }
    }

    pub fn log10_p_error(&mut self, p_log10_error: f64) {
        self.gq((p_log10_error * -10.0).round() as i64)
    }

    pub fn gq(&mut self, gq: i64) {
        self.gq = gq
    }

    pub fn has_likelihoods(&self) -> bool {
        !self.pl.is_empty()
    }

    pub fn pl(&mut self, pl: Vec<i64>) {
        self.pl = pl
    }

    pub fn has_ad(&self) -> bool {
        self.ad.len() > 0
    }

    pub fn has_dp(&self) -> bool {
        self.dp != -1
    }

    pub fn has_gq(&self) -> bool {
        self.gq != -1
    }

    pub fn get_ad(&mut self) -> &mut Vec<i64> {
        &mut self.ad
    }

    pub fn pl_i32(&self) -> Vec<i32> {
        self.pl.iter().map(|i| *i as i32).collect()
    }

    pub fn ad_i32(&self) -> Vec<i32> {
        self.ad.iter().map(|i| *i as i32).collect()
    }

    pub fn no_call_alleles(&mut self, ploidy: usize) {
        self.alleles = vec![ByteArrayAllele::no_call(); ploidy]
    }

    pub fn no_qg(&mut self) {
        self.gq = -1
    }

    pub fn no_annotations(&mut self) {
        self.gq = -1;
        self.ad = Vec::new();
        self.dp = -1;
        self.attributes = HashMap::new();
    }

    pub fn attribute(&mut self, attribute: String, value: AttributeObject) {
        self.attributes.insert(attribute, value);
    }

    pub fn has_attribute(&self, attribute: &String) -> bool {
        self.attributes.contains_key(attribute)
    }

    pub fn get_attribute(&self, attribute: &String) -> Option<&AttributeObject> {
        self.attributes.get(attribute)
    }

    pub fn get_attribute_mut(
        &mut self,
        attribute: String,
        default: AttributeObject,
    ) -> &mut AttributeObject {
        self.attributes.entry(attribute).or_insert(default)
    }

    pub fn alleles(&mut self, alleles: Vec<ByteArrayAllele>) {
        self.alleles = alleles
    }
    // pub fn genotype_likelihood_calculator(&self,)

    // fn calculate_genotype_counts_using_tables_and_validate()

    pub fn get_type(&mut self) -> &GenotypeType {
        if self.genotype_type.is_none() {
            self.genotype_type = Some(self.determine_type())
        }

        return self.genotype_type.as_ref().unwrap();
    }

    /**
     * @return true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this method will return false.
     */
    pub fn is_hom(&mut self) -> bool {
        return self.is_hom_ref() || self.is_hom_var();
    }

    /**
     * @return true if all observed alleles are ref; if any alleles are no-calls, this method will return false.
     */
    pub fn is_hom_ref(&mut self) -> bool {
        return self.get_type() == &GenotypeType::HomRef;
    }

    /**
     * @return true if all observed alleles are alt; if any alleles are no-calls, this method will return false.
     */
    pub fn is_hom_var(&mut self) -> bool {
        return self.get_type() == &GenotypeType::HomVar;
    }

    /**
     * @return true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
     */
    pub fn is_het(&mut self) -> bool {
        return self.get_type() == &GenotypeType::Het;
    }

    /**
     * @return true if we're het (observed alleles differ) and neither allele is reference; if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
     */
    pub fn is_het_non_ref(&mut self) -> bool {
        return (self.get_type() == &GenotypeType::Het)
            && !self.alleles[0].is_ref
            && !self.alleles[1].is_ref;
    }

    /**
     * Internal code to determine the type of the genotype from the alleles vector
     * @return the type
     */
    fn determine_type(&mut self) -> GenotypeType {
        // TODO -- this code is slow and could be optimized for the diploid case
        if self.alleles.is_empty() {
            return GenotypeType::Unavailable;
        }

        let mut saw_no_call = false;
        let mut saw_multiple_alleles = false;
        let mut first_call_allele = None;

        for allele in self.alleles.iter() {
            if allele.is_no_call {
                saw_no_call = true;
            } else if first_call_allele.is_none() {
                first_call_allele = Some(allele);
            } else if allele != first_call_allele.unwrap() {
                saw_multiple_alleles = true;
            }
        }

        if saw_no_call {
            if first_call_allele.is_none() {
                return GenotypeType::NoCall;
            } else {
                return GenotypeType::Mixed;
            }
        };

        match first_call_allele {
            None => {
                panic!("BUG: There are no alleles present in this genotype but the alleles list is not empty");
            }
            Some(first_call_allele) => {
                if saw_multiple_alleles {
                    return GenotypeType::Het;
                } else if first_call_allele.is_ref {
                    return GenotypeType::HomRef;
                } else {
                    return GenotypeType::HomVar;
                }
            }
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct GenotypesContext {
    // sample_names_in_order: Vec<String>,
    genotypes: Vec<Genotype>,
    max_ploidy: i32,
}

impl GenotypesContext {
    pub fn empty() -> GenotypesContext {
        GenotypesContext {
            genotypes: Vec::new(),
            max_ploidy: -1,
        }
    }

    pub fn create(size: usize) -> GenotypesContext {
        GenotypesContext {
            genotypes: Vec::with_capacity(size),
            max_ploidy: -1,
        }
    }

    pub fn new(genotypes: Vec<Genotype>) -> GenotypesContext {
        GenotypesContext {
            // sample_names_in_order: Vec::new(),
            genotypes,
            max_ploidy: -1,
        }
    }

    pub fn contains_sample(&self, sample_name: &String) -> bool {
        return self
            .genotypes
            .par_iter()
            .any(|g| &g.sample_name == sample_name);
    }

    pub fn add(&mut self, genotype: Genotype) {
        self.genotypes.push(genotype)
    }

    pub fn is_empty(&self) -> bool {
        self.genotypes.len() == 0
    }

    pub fn size(&self) -> usize {
        self.genotypes.len()
    }

    pub fn genotypes(&self) -> &Vec<Genotype> {
        &self.genotypes
    }

    pub fn genotypes_mut(&mut self) -> &mut Vec<Genotype> {
        &mut self.genotypes
    }

    pub fn get(&self, index: usize) -> Genotype {
        self.genotypes[index].clone()
    }

    pub fn get_dp(&self) -> i64 {
        self.genotypes[0].dp
    }

    pub fn get_max_ploidy(&mut self, default_ploidy: usize) -> i32 {
        if self.max_ploidy == -1 {
            self.max_ploidy = 0;
            for g in &self.genotypes {
                self.max_ploidy = std::cmp::max(self.max_ploidy, g.ploidy as i32)
            }

            if self.max_ploidy == 0 {
                self.max_ploidy = default_ploidy as i32
            }
        }
        return self.max_ploidy;
    }

    pub fn len(&self) -> usize {
        self.genotypes.len()
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum GenotypeType {
    /** The sample is no-called (all alleles are NO_CALL */
    NoCall,
    /** The sample is homozygous reference */
    HomRef,
    /** The sample is heterozygous, with at least one ref and at least one one alt in any order */
    Het,
    /** All alleles are non-reference */
    HomVar,
    /** There is no allele data availble for this sample (alleles.isEmpty) */
    Unavailable,
    /** Some chromosomes are NO_CALL and others are called */
    Mixed,
}
