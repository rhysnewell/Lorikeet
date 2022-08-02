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
      fields=['variants/CHROM', 'variants/POS', 'variants/DP', 'calldata/GT', 'calldata/AD', 'calldata/DP'])

vcf['calldata/AD'][0, 0, :]
haplotype_arrays = []
for j in range(vcf['calldata/AD'].shape[1]): # samples (or populations)
    population_j = []
    max_depth = vcf['calldata/DP'][:, j].max() # All variants have to info for this many individuals
    for i in range(vcf['calldata/AD'].shape[0]): # variant sites
        site_i = np.array([-1 for _ in range(max_depth)])
        depth_count = 0
        for k in range(vcf['calldata/AD'].shape[2]): # alleles at site in sample
            if k <= 0:
                continue

            ad = vcf['calldata/AD'][i, j, k] # allele depth
            site_i[depth_count:(depth_count + ad)] = k
            depth_count += ad # update minimum index

        population_j.append(site_i)
    haplotype_arrays.append(allel.HaplotypeArray(population_j))



col_dict = ['sample_1', 'sample_2']
for variant_i in range(vcf['calldata/AD'].shape[0]):
    col_dict.append(f'variant_{variant_i}')

fst_df = np.zeros(shape=(int(scipy.special.binom(vcf['calldata/DP'].shape[1], 2)) * 2, vcf['calldata/AD'].shape[0] + 2))
mean_fst_df = np.zeros(shape=(vcf['calldata/DP'].shape[1], vcf['calldata/DP'].shape[1]))
row = 0
for sample1 in range(vcf['calldata/AD'].shape[1]):
    for sample2 in range(vcf['calldata/AD'].shape[1]):
        if sample1 == sample2:
            continue

        num, den = allel.hudson_fst(haplotype_arrays[sample1], haplotype_arrays[sample2])
        fst = (num / den)
        np.nan_to_num(fst, copy=False, nan=0.0, posinf=None, neginf=None)
        fst[fst < 0] = 0

        mean_fst = max(fst.mean(), 0)
        mean_fst_df[sample1, sample2] = mean_fst
        fst_df[row, 0:2] = [sample1, sample2]
        fst_df[row, 2:] = fst
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