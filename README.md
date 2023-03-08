# Sperm transmission mosaicism

This repository collects pipeline, code, and some intermediate results for the study of mosaic SNV/Indels for sperm, blood, and saliva samples of a small cohort. Raw WGS data of this study is available on SRA under [PRJNA753973](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA753973/). Genotyping results of each sample is provided [here](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Analysis/Genotyping_in_different_samples.csv). Re-analysis of the transmission of paternal mosaic variants from the ASD trios is provided [here](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Analysis/Genotyping_ASD_trios.csv).

The pipelines and analysis are derived from our recent [large-scale sperm study](https://www.sciencedirect.com/science/article/pii/S0092867421008837), of which the code could be found [here](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism).

-----------------------------------

### 1. Pipelines for the process of whole-genome sequencing data

#### 1.1 Pipelines for WGS data process and quality control

[Pipelines](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Pipelines/Preprocessing) for pre-processing of the bams.

Code for [depth of coverage](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Depth_of_coverage.r) and [insertsize distribution](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Insert_size.r).

#### 1.2 Code for the population origin analysis

[Pipeline](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Whole_genome_variant_extraction_and_PCA.sh) for population analysis, and [code](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Plot_PCA.r) for plot.

#### 1.3 Pipelines for mosaic SNV/indel calling and variant annotations

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode) for MuTect2 (single mode) with the "Full Panel of Normal" version is used. The MuTect2 (single mode) result is followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

-----------------------------------

### 2. Pipelines for the process of Massive Parallel Amplicon Sequencing (MPAS)
#### 2.1 Pipelines for MPAS data alignment and processing
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/MPAS_and_snMPAS_processing_pipeline) for alignment, processing, and germline variant calling of MPAS reads.

#### 2.2 Pipelines for AF quantification and variant annotations

[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

[Code](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Analysis/transmission_ampliseq_analysis_plotting.py) to filter and annotate on MPAS data.

-----------------------------------

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting

After variant calling from different strategies, variants were annotated and filtered by [a python script](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Analysis/transmission_ampliseq_analysis_plotting.py) and positive mosaic variants as well as the corresponding transmission to multiple samples and additional information were annotated.

#### 3.2 Pipelines for statistically analysis, and the related plotting

[Code](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Plots/transmission_ampliseq_analysis_plotting.py) for the estimation of expected transmissions assuming independent transmission via a dynamaic programming algorithm.

[Code](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Analysis/Simulation/Plot.r) and [example data](https://github.com/shishenyxx/Sperm_transmission_mosaicism/tree/main/Analysis/Simulation) for the permutation analysis to estimate the indepence of transmission in each family.

[Code](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Plots/transmission_WGS_ASD_trios.py) and [data](https://github.com/shishenyxx/Sperm_transmission_mosaicism/blob/main/Analysis/Genotyping_ASD_trios.csv) for the re-analysis of transmission in 8 ASD families previously analyzed in the [first](https://www.nature.com/articles/s41591-019-0711-0) and [second](https://doi.org/10.1016/j.cell.2021.07.024) study.

-----------------------------------

### 4. Contact:

:email: Martin Breuss: [martin.breuss@cuanschutz.edu](mailto:martin.breuss@cuanschutz.edu)

:email: Xiaoxu Yang: [xiy010@health.ucsd.edu](mailto:xiy010@health.ucsd.edu), [yangxiaoxu-shishen@hotmail.com](mailto:yangxiaoxu-shishen@hotmail.com)

:email: Joseph Gleeson: [jogleeson@health.ucsd.edu](mailto:jogleeson@health.ucsd.edu), or the Gleeson lab [gleesonlab@health.ucsd.edu](gleesonlab@health.ucsd.edu)

-----------------------------------

### 5. Cite the code
Breuss MW, Yang X, <i>et al.</i>, Gleeson JG. [Unbiased mosaic variant assessment in sperm: a cohort study to test predictability of transmission.](https://elifesciences.org/articles/78459) 2022. (<i>eLife</i>, DOI:[10.7554/eLife.78459](https://doi.org/10.7554/elife.78459), PMID:[35787314](https://pubmed.ncbi.nlm.nih.gov/35787314/))

<img src="https://user-images.githubusercontent.com/17311837/223878671-d6e925a7-d834-41ac-9728-5403ec23662a.png" alt="Sperm_Mosaic_Cover" width=80%> 
