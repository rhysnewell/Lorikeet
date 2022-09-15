use pyo3::prelude::*;
use pyo3::types::IntoPyDict;

pub fn calculate_fst(
    output_prefix: &str,
    genome_name: &str,
    vcf_path: &str,
    ploidy: usize,
    depth_per_sample_filter: i64
) -> PyResult<()> {
    pyo3::prepare_freethreaded_python();
    Python::with_gil(|py| {
        // let sys = py.import("allel")?;
        // let version: String = sys.getattr("version")?.extract()?;
        let np = py.import("numpy")?;
        let allel = py.import("allel")?;
        let pr = py.import("polars")?;
        let depth_per_sample_str = format!("{depth_per_sample_filter}");
        let locals = [
            ("output_prefix", output_prefix),
            ("genome_name", genome_name),
            ("vcf_path", vcf_path),
            ("depth_per_sample", depth_per_sample_str.as_str())
        ].into_py_dict(py);
        let code =r#"
import polars as pr
import numpy as np
import allel
import scipy.special
import warnings

# Warnings break the lorikeet interface
warnings.filterwarnings('ignore')
depth_per_sample = int(depth_per_sample)
try:
    vcf = allel.read_vcf(vcf_path,
          fields=['variants/CHROM', 'variants/POS', 'variants/DP', 'variants/QF', 'calldata/GT', 'calldata/AD', 'calldata/DP'])
except FileNotFoundError:
    vcf = allel.read_vcf(f"{vcf_path}.gz",
          fields=['variants/CHROM', 'variants/POS', 'variants/DP', 'variants/QF', 'calldata/GT', 'calldata/AD', 'calldata/DP'])

col_dict = ['sample_1', 'sample_2']
n_variants = 0
for variant_i in range(vcf['calldata/AD'].shape[0]):
    if vcf['variants/QF'][variant_i] == 'false':
            continue # do not include unqualified variants
    col_dict.append(f'variant_{variant_i}')
    n_variants += 1

# fst_df = np.zeros(shape=(int(scipy.special.binom(vcf['calldata/DP'].shape[1], 2)), n_variants + 2))
mean_fst_df = np.zeros(shape=(vcf['calldata/DP'].shape[1], vcf['calldata/DP'].shape[1]))
row = 0
allele_counts = vcf['calldata/AD'][vcf['variants/QF'] == 'true']
population_size = vcf['calldata/DP'][vcf['variants/QF'] == 'true']
# allele_counts = vcf['calldata/AD']
# population_size = vcf['calldata/DP']

for sample1 in range(vcf['calldata/AD'].shape[1]):
    for sample2 in range(sample1 + 1):
        if sample1 == sample2:
            continue

        sample1_depths = population_size[:, sample1]
        sample1_allele_counts = allele_counts[:, sample1, (allele_counts!=-1)[:, sample1].all(axis=0)]
        # sample2_depths = np.rot90(np.repeat([population_size[:, sample2]], 2, axis=0))
        sample2_depths = population_size[:, sample2]
        sample2_allele_counts = allele_counts[:, sample2, (allele_counts!=-1)[:, sample2].all(axis=0)]

        variants_to_include_1 = (sample1_depths >= depth_per_sample) #+ ((sample1_allele_counts / sample1_depths) > 0.01)
        variants_to_include_2 = (sample2_depths >= depth_per_sample) #+ ((sample2_allele_counts / sample2_depths) > 0.01)
        variants_to_include = np.array([variants_to_include_1, variants_to_include_2]).all(axis=0)

        h1 = allel.HaplotypeArray(sample1_allele_counts[variants_to_include])
        h2 = allel.HaplotypeArray(sample2_allele_counts[variants_to_include])

        num, den = allel.hudson_fst(h1, h2)
        fst = (num / den)

        np.nan_to_num(fst, copy=False, nan=0.0, posinf=None, neginf=None)

        fst[fst < 0] = 0.0
        fst[fst > 1] = 1.0
        #fst_df[row, 2:] = fst

        #mean_fst = max(fst_df[row, 2:].mean(), 0)
        mean_fst_df[sample1, sample2] = fst.mean()
        mean_fst_df[sample2, sample1] = fst.mean()
        #fst_df[row, 0:2] = [sample1, sample2]
        row += 1

#pr_df = pr.DataFrame(fst_df, columns=col_dict)
samples = [str(i + 1) for i in range(vcf['calldata/DP'].shape[1])]
np.nan_to_num(mean_fst_df, copy=False, nan=0, posinf=0, neginf=0)
mean_fst_df = np.insert(mean_fst_df, 0, np.array([i + 1 for i in range(vcf['calldata/DP'].shape[1])]), axis=1)
samples.insert(0, "SampleID")
mean_fst_df = pr.DataFrame(mean_fst_df, columns=samples)
#pr_df.write_csv(file=f"{output_prefix}/{genome_name}_fst_values.tsv", sep='\t')
mean_fst_df.write_csv(file=f"{output_prefix}/{genome_name}_sample_fst_values.tsv", sep='\t')"#;
        py.run(code, None, Some(&locals))?;
        Ok(())
    })
}