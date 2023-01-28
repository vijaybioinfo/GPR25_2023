# GPR25_2023


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
> Relevant scripts all located in: ./figures


Downstream Analysis
------------
* DGEA - You can follow [our in-house pipeline](https://github.com/vijaybioinfo/dgea)
> Relevant scripts all located in: ./downstream_analysis


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
