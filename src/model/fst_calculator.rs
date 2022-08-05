use pyo3::prelude::*;
use pyo3::types::IntoPyDict;

pub fn calculate_fst(output_prefix: &str, genome_name: &str, ploidy: usize) -> PyResult<()> {
    pyo3::prepare_freethreaded_python();
    Python::with_gil(|py| {
        // let sys = py.import("allel")?;
        // let version: String = sys.getattr("version")?.extract()?;
        let np = py.import("numpy")?;
        let allel = py.import("allel")?;
        let pr = py.import("polars")?;
        let locals = [
            ("output_prefix", output_prefix),
            ("genome_name", genome_name)
        ].into_py_dict(py);
        let code =r#"
import polars as pr
import numpy as np
import allel
import scipy.special

vcf = allel.read_vcf(f'{output_prefix}/{genome_name}.vcf',
      fields=['variants/CHROM', 'variants/POS', 'variants/DP', 'variants/QF', 'calldata/GT', 'calldata/AD', 'calldata/DP'])

col_dict = ['sample_1', 'sample_2']
n_variants = 0
for variant_i in range(vcf['calldata/AD'].shape[0]):
    if vcf['variants/QF'][variant_i] == 'false':
            continue # do not include unqualified variants
    col_dict.append(f'variant_{variant_i}')
    n_variants += 1

fst_df = np.zeros(shape=(int(scipy.special.binom(vcf['calldata/DP'].shape[1], 2)) * 2, n_variants + 2))
mean_fst_df = np.zeros(shape=(vcf['calldata/DP'].shape[1], vcf['calldata/DP'].shape[1]))
row = 0
allele_counts = vcf['calldata/AD'][vcf['variants/QF'] == 'true']

for sample1 in range(vcf['calldata/AD'].shape[1]):
    for sample2 in range(vcf['calldata/AD'].shape[1]):
        if sample1 == sample2:
            continue
        for var_i in range(allele_counts.shape[0]):
            h1 = allel.HaplotypeArray([allele_counts[var_i, sample1][allele_counts[var_i, sample1] != -1]])
            h2 = allel.HaplotypeArray([allele_counts[var_i, sample2][allele_counts[var_i, sample2] != -1]])
            num, den = allel.hudson_fst(h1, h2)
            fst = (num / den)

            np.nan_to_num(fst, copy=False, nan=0.0, posinf=None, neginf=None)

            fst[fst < 0] = 0
            fst_df[row, var_i + 2] = fst

        mean_fst = max(fst_df[row, 2:].mean(), 0)
        mean_fst_df[sample1, sample2] = mean_fst
        fst_df[row, 0:2] = [sample1, sample2]
        row += 1

pr_df = pr.DataFrame(fst_df, columns=col_dict)
samples = [str(i + 1) for i in range(vcf['calldata/DP'].shape[1])]
mean_fst_df = np.insert(mean_fst_df, 0, np.array([i + 1 for i in range(vcf['calldata/DP'].shape[1])]), axis=1)
samples.insert(0, "SampleID")
mean_fst_df = pr.DataFrame(mean_fst_df, columns=samples)
pr_df.write_csv(file=f"{output_prefix}/{genome_name}_fst_values.tsv", sep='\t')
mean_fst_df.write_csv(file=f"{output_prefix}/{genome_name}_sample_fst_values.tsv", sep='\t')"#;
        py.run(code, None, Some(&locals))?;
        Ok(())
    })
}