# GPR25_2023

------------

Description
------------

Tissue-resident memory (TRM) CD8+ T cells are key players that orchestrate protective anti-viral and anti-tumor immune responses. The molecules that support their development and function are not fully defined. Here, we report on the regulation and function of an orphan G-protein coupled receptor, GPR25 that is expressed at high levels in TRM cells compared to non-TRM cells. TGF-B, a key cytokine involved in the development of TRM cells, induces the expression of GPR25 in CD8+ T cells, and TRM-associated cis-regulatory elements in the GPR25 locus show binding of SMAD1, a key transcription factor downstream of TGF-B signaling. Using Gpr25-deficient T cells in an LCMV infection model, we show that Gpr25 acts in a cell-intrinsic manner to promote the development of both primary and secondary TRM cells in the liver. Gpr25 deficiency in T cells also impairs their capacity to develop into lung TRM cells and control lung metastases. Single-cell transcriptomic analysis of Gpr25-deficient memory T cells and TRM cells showed potential defects in restraining the ZEB2-S1PR5 pathway that drives tissue egress and in pathways supporting a stem-like memory program as opposed to effector differentiation. Our findings support the concept that modulating Gpr25 function may be an attractive therapeutic option to boost the magnitude and quality of TRM responses generated in the context of infection and cancer. 

This repository contains the data and scripts used to analyze the aforementioned samples.

Requirements
------------

This project was done using the following modules/programs:

* [R](https://cran.r-project.org/) (v4.0.1)
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v3.1.0)
* [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)
* [Seurat](https://satijalab.org/seurat) (v4.1.1)


Raw data
------------
The single-cell RNA-seq raw and processed files can be downloaded through the following GEO accession number: [GSE223627](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223627). 

Raw data pre-processing  
------------

### Single-cell
* To do the 10x demultiplexing and mapping just pull [our in-house pipeline](https://github.com/vijaybioinfo/cellranger_wrappeR) using Cell Ranger.
* To do the donor the multiplexing just pull [our in-house pipeline](https://github.com/vijaybioinfo/ab_capture).
* To do the single-cell quality control just pull [our in-house pipeline](https://github.com/vijaybioinfo/quality_control).
* To do the doublet detection use [our in-house pipeline](https://github.com/vijaybioinfo/doublet_detection) using Scrublet. 
* To generate the clustering of single-cell data just pull [our in-house pipeline](https://github.com/vijaybioinfo/clustering) using Seurat.

> Relevant scripts all located in: ./pre-processing  

For more specific information about the data generation and processing, please check the "methods" section within the manuscript.  


Figures
------------
> Relevant scripts all located in: ./figures_all.R

Downstream Analysis
------------
* DGEA - You can follow [our in-house pipeline](https://github.com/vijaybioinfo/dgea)
> Relevant script all located in: ./downstream_liver

Citation
--------------

Please cite the following manuscript if you are using this repository:


Maintainers
-----------

Current maintainers:
* Francisco Emmanuel Castaneda Castro (fcastaneda@lji.org) 
* Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org)

Vijayanand Lab.  
Division of Vaccine Discovery La Jolla Institute for Immunology La Jolla, CA 92037, USA


Contact
-----------
Please email Francisco Emmanuel Castaneda Castro (fcastaneda@lji.org), Ciro Ramírez-Suástegui (ksuasteguic@gmail.com, ciro@lji.org) and/or Vijayanand Pandurangan (vijay@lji.org).
