use ordered_float::{NotNan, OrderedFloat};
use ndarray::{Array, Array2, ArrayBase, OwnedRepr};
use model::genotype_allele_counts::GenotypeAlleleCounts;
use model::genotype_likelihood_calculator::GenotypeLikelihoodCalculator;
use model::variants::Allele;

pub enum GenotypeAssignmentMethod {
    BestMatchToOriginal,
    DoNotAssignGenotypes,
    SetToNoCall,
    SetToNoCallNoAnnotations,
    UsePLsToAssign,
}

pub struct Genotype {
    pub ploidy: usize,
    pub likelihoods: Vec<OrderedFloat<f64>>,
    pub alleles: Vec<Allele>,
}

impl Genotype {
    pub fn build(default_ploidy: usize, likelihoods: Vec<OrderedFloat<f64>>) -> Genotype {
        Genotype {
            ploidy: default_ploidy,
            alleles: Vec::with_capacity(likelihoods.len()),
            likelihoods,
        }
    }

    pub fn build_from_alleles(alleles: Vec<Allele>) -> Genotype {
        Genotype {
            ploidy: alleles.len(),
            likelihoods: vec![OrderedFloat(0.0); alleles.len()],
            alleles,
        }
    }

    pub fn get_ploidy(&self) -> usize { self.ploidy }

    pub fn get_likelihoods(&mut self) -> &mut Vec<OrderedFloat<f64>> {
        &mut self.likelihoods
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

    pub fn genotypes(self) -> Vec<Genotype> {
        self.genotypes
    }

    pub fn get(&self, index: usize) -> Genotype {
        self.genotypes[index].clone()
    }
}