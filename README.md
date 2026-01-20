## Optimizing Global Genomic Surveillance for Early Detection of Emerging SARS-CoV-2 Variants

<!-- *Haogao Gu#, Jifan Li#, Wanying Sun, Mengting Li, Kathy Leung, Joseph T. Wu, Hsiang-Yu Yuan, Matthew R. McKay, Ning Ning, Leo L.M. Poon\** -->

### Project Overview:

The global spread of viruses highlights the need for timely and effective genomic surveillance to detect new variants and inform rapid public health responses. However, high costs and uneven sequencing capacity hinder equitable global implementation. Surveillance focused on international travelers at major travel hubs has been proposed as a way to complement robust local surveillance, but its potential benefits have not been fully quantified. Here, we develop and calibrate a multiple-strain metapopulation model of global SARS-CoV-2 transmission using extensive epidemiological, phylogenetic, and high-resolution air travel data. Retrospective analyses of the Omicron BA.1/BA.2 emergence and forward simulations for hypothetical novel variants show that targeted enhancement of traveler surveillance at key hubs can shorten variant detection delays, with reduced total surveillance efforts. Practical “non-disruptive” strategies, such as prioritizing a small number of highly connected hubs, consistently outperformed baseline approaches and remained effective across a range of variant transmissibility and vaccine effectiveness scenarios. These results provide a quantitative framework for strengthening global genomic surveillance through targeted, complementary strategies that preserve local capacity while improving preparedness for future pandemics.

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
