# Code_for_Wang_et_al2020
This folder contains R script and RData to generate the main Figures from Wang et al. (2020): Theoretical and empirical quantification of the accuracy of polygenic scores in ancestry divergent populations

Source code of a custom C++ used to calculate LD-based statistics
Compile C++ code
> make

> ./ldcorpair --help

Options
> --bfile      : Binary PLINK format for genotypes.

> --snplist    : Specify the list of SNPs containing the target SNPs.

> --window-kb  : Specify the window size to calculate LD correlations. Default is 100 kb.

> --out        : A prefix for the output file [prefix].ldcor.pairs.

[Note] Missing values are imputed to the major allele.

Overview of output

> CHR	SNP1	SNP2	POS1_KB	POS2_KB	A1_1	A1_2	A2_1	A2_2	FREQ_A1_1	FREQ_A1_2	R_LD	R_LD_SQ	R_LD_SQ_CORR

> 1	rs28536514	rs58686784	910.409	810.78	C	C	T	G	0.19062	0.304841	0.103561	0.0107248	0.00920504

> 1	rs28536514	rs114224245	910.409	811.136	C	C	T	G	0.19062	0.0234493	0.00218953	4.79404e-06	-0.00151496


The annotated script: "Annotated_script_to_calculate_expected_relative_acccuracy.R" was added in Dec 2021 and used in Okbay et al. (2021) to calculate the expected relative accuracy of a Polygenic Score predicting educational attainment in a European ancestries sample relative to an African ancestries sample. The same method was used in Wang et al. (2020).
