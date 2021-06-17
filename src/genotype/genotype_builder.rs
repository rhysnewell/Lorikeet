use genotype::genotype_allele_counts::GenotypeAlleleCounts;
use genotype::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use model::variants::Allele;
use genotype::genotype_likelihoods::GenotypeLikelihoods;
use std::collections::HashMap;

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

#[derive(Debug, Clone)]
pub struct Genotype {
    pub ploidy: usize,
    pub pl: GenotypeLikelihoods,
    pub alleles: Vec<Allele>,
    pub ad: Vec<i64>,
    pub dp: i64,
    pub gq: i64,
    pub is_phased: bool,
    pub attributes: HashMap<String, Vec<f64>>
}

impl Genotype {
    const HAPLOID_NO_CALL: Vec<Allele> = vec![Allele::fake(false)];
    const DIPLOID_NO_CALL: Vec<Allele> = vec![Allele::fake(false); 2];


    pub fn build(default_ploidy: usize, likelihoods: Vec<f64>) -> Genotype {
        Genotype {
            ploidy: default_ploidy,
            alleles: Vec::with_capacity(likelihoods.len()),
            pl: GenotypeLikelihoods::from_log10_likelihoods(likelihoods),
            dp: -1,
            gq: -1,
            ad: Vec::with_capacity(likelihoods.len()),
            is_phased: false,
            attributes: HashMap::new(),
        }
    }

    pub fn build_from_alleles(alleles: Vec<Allele>) -> Genotype {
        Genotype {
            ploidy: alleles.len(),
            pl: GenotypeLikelihoods::from_log10_likelihoods(vec![0.0; alleles.len()]),
            dp: -1,
            gq: -1,
            ad: Vec::with_capacity(alleles.len()),
            is_phased: false,
            attributes: HashMap::new(),
            alleles,
        }
    }

    pub fn get_ploidy(&self) -> usize { self.ploidy }

    pub fn get_likelihoods(&mut self) -> &mut Vec<f64> {
        &mut self.pl.get_as_vector()
    }

    pub fn num_likelihoods(&mut self, num_alleles: i64, ploidy: i64) -> i64 {
        self.pl.num_likelihoods(num_alleles, ploidy)
    }

    pub fn log10_p_error(&mut self, p_log10_error: f64) {
        self.gq((p_log10_error * -10.0) as i64)
    }

    pub fn gq(&mut self, gq: i64) {
        self.gq = gq
    }

    pub fn has_likelihoods(&self) -> bool {
        !self.pl.is_empty()
    }

    pub fn pl(&mut self, pl: GenotypeLikelihoods) {
        self.pl = pl
    }

    pub fn has_ad(&self) -> bool {
        self.ad.len() == self.alleles.len()
    }

    pub fn get_ad(&mut self) -> &mut Vec<i64> {
        &mut self.ad
    }

    pub fn no_call_alleles(&mut self, ploidy: usize) {
        self.alleles = vec![Allele::no_call(); ploidy]
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

    pub fn attribute(&mut self, attribute: String, value: Vec<f64>) {
        self.attributes.insert(attribute, value);
    }

    pub fn has_attribute(&self, attribute: &String) -> bool {
        self.attributes.contains_key(attribute)
    }

    pub fn get_attribute(&mut self, attribute: &String) -> &mut Vec<f64> {
        self.attributes.entry(attribute.clone()).or_insert(vec![std::f64::NAN; self.alleles.len()])
    }

    // pub fn genotype_likelihood_calculator(&self,)


    // fn calculate_genotype_counts_using_tables_and_validate()
}

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

    pub fn get(&self, index: usize) -> Genotype {
        self.genotypes[index].clone()
    }

    pub fn get_dp(&self) -> i64 {
        self.genotypes[0].dp
    }

    pub fn get_max_ploidy(&mut self, default_ploidy: usize) -> i32 {
        if self.max_ploidy == -1 {
            self.max_ploidy = 0;
            for g in self.genotypes {
                self.max_ploidy = std::cmp::max(self.max_ploidy, g.ploidy as i32)
            }

            if self.max_ploidy == 0 {
                self.max_ploidy = default_ploidy as i32
            }
        }
        return self.max_ploidy
    }
}