# srac2lrac
R code for calculating the concordance statistic between short read derived MAG and long read chromosomes, as described in preprint entitled "Analysis procedures for assessing recovery of high quality, complete, closed genomes from Nanopore long read metagenome sequencing" (Arumugam, Bessarab, Haryono, Liu, Zuniga-Montanez, Roy, Qiu, Drautz-Moses, Law, Wuertz, Lauro, Huson and Williams, 2020)

The following files are provided:
1. srac2lrac.R: a working set of R functions (not yet in package form) for computing the concordance statistic. The input is a data frame containing the tabular output of a BLASTN analyis of LRAC sequences (queries) against SRAC sequences (subjects).
2. local-plot-functions.R: a set of functions for visualising the output. The function "summary.plot.for.srac2lrac.single.ext3" was used to generate Figure 1 and Supplementary Figures 2-34.
3. An analysis from the data set denoted PAO2 in the above paper. Here "srac_spades2corrlrac_canu.out" is the tabular output from a BLASTN analysis in which SRAC sequences (contigs from short read assembly) are treated as queries and LRAC sequences (contigs from long read assembly) are treated as subjects. This data should be read into R as a data frame.
4. There are various other bits of data needed. A data frame that provides summary data for the assembly, indexed by contig identifier, with the following columns: cov [coverage], len [length, in bp] and gc [GC content of contig]. For SPAdes assemblies we extract this data from the assembly contig file using our R package RKXM. The RKXM function is used to generate a reduced version of this table (including a cutoff for minimal contig length, and converting . 1) .
4. The primary output from the srac2lrac.r code is from the function , which returns a table ('blastn7' in the above example) with the following columns: . Using this table you can assess pairs of LRAC and short read bins are likely to be putative congnates. The example plots used in the paper (Figure 1 and Supplementary Figures 2-34 let you explore the support for the concordance statistic results in more detail.
