# srac2lrac
R code for calculating the concordance statistic between short read derived MAG and long read chromosomes, as described in preprint entitled "Analysis procedures for assessing recovery of high quality, complete, closed genomes from Nanopore long read metagenome sequencing" (Arumugam, Bessarab, Haryono, Liu, Zuniga-Montanez, Roy, Qiu, Drautz-Moses, Law, Wuertz, Lauro, Huson and Williams, 2020)

The following files are provided:
1. srac2lrac.R: a set of R functions (not yet in package form) for computing the concordance statistic. The input is a data frame containing the tabular output of a BLASTN analyis of LRAC sequences (queries) against SRAC sequences (subjects).
2. local-plot-functions.R: a set of functions for visualising the output. The function "summary.plot.for.srac2lrac.single.ext3" was used to generate Figure 1 and Supplementary Figures 2-34.
3. Analysis logs from each of the four data sets analysed in the paper. These are provided as-is, and filepaths will have to be modified according to your local configuration. For continuity all output filenames are the same as those provided in our associated Zenodo submission.
