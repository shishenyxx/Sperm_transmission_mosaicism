# Sperm_transmission_mosaicism

This repository collects pipelines, codes, and some intermediate results for the study of mosaic SNV/Indels for sperm, blood, and saliva samples of a small cohort. Raw data of this study is available on SRA under []().

-----------------------------------

### 1. Pipelines for the process of whole-genome sequencing data

#### 1.1 Pipelines for WGS data process and quality control

[Pipelines](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Pipelines/Preprocessing) for pre-processing of the bams.

Codes for [depth of coverage](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Depth_of_coverage.r) and [insertsize distribution](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Insert_size.r).

#### 1.2 Codes for the population origin analysis

[Pipeline](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Whole_genome_variant_extraction_and_PCA.sh) for population analysis, and [codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Plot_PCA.r) for plot.

#### 1.3 Pipelines for mosaic SNV/indel calling and variant annotations

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode) for MuTect2 (single mode) with the "Full Panel of Normal" version is used. The MuTect2 (single mode) result is followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

-----------------------------------

### 2. Pipelines for the process of Massive Parallel Amplicon Sequencing (MPAS)
#### 2.1 Pipelines for MPAS data alignment and processing
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/MPAS_and_snMPAS_processing_pipeline) for alignment, processing, and germline variant calling of MPAS reads.

#### 2.2 Pipelines for AF quantification and variant annotations

[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

[Codes]() to filter and annotate on MPAS data.

-----------------------------------

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting

After variant calling from different strategies, variants were annotated and filtered by [a python script]() and positive mosaic variants as well as the corresponding transmission to multiple samples and additional information were annotated.

#### 3.2 Pipelines for statistically analysis, and the related plotting

[Codes]() for the estimation of expected transmissions assuming independent transmission via a dynamaic programming algorithm.


