## Optimizing Global Genomic Surveillance for Early Detection of Emerging SARS-CoV-2 Variants

<!-- *Haogao Gu#, Jifan Li#, Wanying Sun, Mengting Li, Kathy Leung, Joseph T. Wu, Hsiang-Yu Yuan, Matthew R. McKay, Ning Ning, Leo L.M. Poon\** -->

### Project Overview:

The global spread of viruses presents significant challenges to public health, underscoring the need for effective monitoring and control measures. Genomic surveillance, which involves the sequencing and analysis of viral genomes, is crucial for detecting new variants and guiding interventions. However, high costs and uneven resource distribution hinder its global implementation.

This project demonstrates that optimizing genomic surveillance by focusing on high-risk groups—specifically international travelers and major travel hubs—can significantly enhance the early detection of emerging SARS-CoV-2 variants. Utilizing a comprehensive metapopulation multiple-strain model calibrated with extensive data, we show that reallocating existing resources towards these key areas reduces detection times without additional costs. This approach offers a sustainable solution for global surveillance networks, providing actionable insights for policymakers to strengthen global health security and improve preparedness for future pandemics.

### Instructions for use:

Readme files are available in each subfolder to guide the workflow. The project is organized into the following subfolders:
- [data](data/): the raw data used in the project. Note that some data are not publicly available and readers should refer to their original source for access.
- [scripts](scripts/): the scripts used for data processing, model building, model fitting, and model simulation.
- [results](results/): the results generated from the analysis.
- Reproduction instructions are provided in the following sections. Also refer to the in-line comments in the scripts for detailed instructions.

### System requirements

- R environment:
  - R version 
	```r
	R version 4.4.1 (2024-06-14)
	Platform: x86_64-pc-linux-gnu
	Running under: Ubuntu 22.04.4 LTS

	Matrix products: default
	BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
	LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

	locale:
	[1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
	[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

	time zone: Asia/Hong_Kong
	tzcode source: system (glibc)

	attached base packages:
	[1] stats     graphics  grDevices utils     datasets  methods   base     

	loaded via a namespace (and not attached):
	[1] compiler_4.4.1 cli_3.6.3      jsonlite_1.8.8 rlang_1.1.4  
	```
  - Required R packages:
	please refer to [scripts/data_processing/install_prerequisite.R](scripts/data_processing/install_prerequisite.R) for the list of required R packages.

### Installation and running time

- Model calibration and simulation are the most time-consuming steps in the project. They were run on a high-performance computing cluster with >30 computing nodes, each with 64 threads. The total running time for the model calibration step is approximately 50 days (including debugging and optimization), and the total running time for the model simulation step is approximately 10 days. The running time may vary depending on the computing resources available.
- Installation and compilation of the model code takes relatively less time. it can be done on a personal computer within 1 hour.

### Demo

- Essential intermediate results are provided in the [results/model_data](results/model_data/).
- To reproduce the results and figures presented in the manuscript, please refer to the [analyzing scripts](scripts/final_analysis) and directly run the codes by loading the above intermediate results.
- Expected output should be exactly the same as the results presented in the manuscript.
- Running time for the analyzing scripts is approximately 1-3 hour (most time-consuming part is the bootstrapping step) on a normal personal computer.

### Citation
Cite as:	
- arXiv:2502.00934 [q-bio.PE] (or arXiv:2502.00934v2 [q-bio.PE] for this version)
- https://doi.org/10.48550/arXiv.2502.00934
- ```bibtex
	@misc{gu2025optimizingglobalgenomicsurveillance,
		title={Optimizing Global Genomic Surveillance for Early Detection of Emerging SARS-CoV-2 Variants}, 
		author={Haogao Gu and Jifan Li and Wanying Sun and Mengting Li and Kathy Leung and Joseph T. Wu and Hsiang-Yu Yuan and Maggie H. Wang and Bingyi Yang and Matthew R. McKay and Ning Ning and Leo L. M. Poon},
		year={2025},
		eprint={2502.00934},
		archivePrefix={arXiv},
		primaryClass={q-bio.PE},
		url={https://arxiv.org/abs/2502.00934}, 
	}
	```


---
