# srac2lrac
R code for calculating the concordance statistic between short read derived MAG and long read chromosomes, as described in preprint entitled "Analysis procedures for assessing recovery of high quality, complete, closed genomes from Nanopore long read metagenome sequencing" (Arumugam, Bessarab, Haryono, Liu, Zuniga-Montanez, Roy, Qiu, Drautz-Moses, Law, Wuertz, Lauro, Huson and Williams, 2020)

The following files are provided:
1. srac2lrac.R: a working set of R functions (not yet in package form) for computing the concordance statistic. The input is a data frame containing the tabular output of a BLASTN analyis of LRAC sequences (queries) against SRAC sequences (subjects).
2. local-plot-functions.R: a set of functions for visualising the output. The function "summary.plot.for.srac2lrac.single.ext3" was used to generate Figure 1 and Supplementary Figures 2-34.
3. An analysis from the data set denoted PAO2 in the above paper (see 'srac2lrac-example.R'). Here "srac_spades2corrlrac_canu.out" is the tabular output from a BLASTN analysis in which SRAC sequences (contigs from short read assembly) are treated as queries and LRAC sequences (contigs from long read assembly) are treated as subjects. This data should be read into R as a data frame.
4. There are various other bits of data needed. A data frame that provides summary data for the assembly, indexed by contig identifier, with the following columns: cov [coverage], len [length, in bp] and gc [GC content of contig]. For SPAdes assemblies we extract this data from the assembly contig file using our R package RKXM. The RKXM function is used to generate a reduced version of this table (including a cutoff for minimal contig length, and converting coverage to log10(coverage). An R (bin identifier) named list is required, whose elements are character vectors containing the memberships (by contig identifiers) of a bin. In the example data, this is provided as 'pao2.srBin.rda' (generated directly from the Metabat2 output using the RKXM function 'make.FASTA.file.set'.
5. The primary output from the srac2lrac.r code is from the function 'augment.kappa.table', which returns a table ('blastn7' in the above example) with the following columns: 1) tag: a string formed by <lrac.id>-<sracBin.id>; 2) lrac: the LRAC identifier; 3) sracBin: the identifier for a short read bin; 4) kappa: the kappa statistics; the last four columns specifiy the component statistics of kappa, namely 5) num.srac (p_sarc); 6) pident (\widehat{pid}); 7) al2ql (\widehat{al2ql}) and 8) p_aln. Using this table you can assess pairs of LRAC and short read bins are likely to be putative congnates. The example plots used in the paper (Figure 1 and Supplementary Figures 2-34 let you explore the support for the concordance statistic results in more detail.
6. The file called 'local-plot-functions.r' contains a set of R functions for plotting. This file contains many 'experiments' but the current approach that we use in the preprint is from the function called 'summary.plot.for.srac2lrac.single.ext3'. We have provided examples that produce Supplementary Figure 12 and 13 from the preprint (see 2 pdf files). For Panel C of these plots you will need a reduced form of the per-base coverage data, obtained from mapping both LR and SR data to the set of genomes. Fro the data in the preprint, we have provided these in the file called 'reducedCoverageProfiles.RData', along with the R code that genrated them in 'cov-log-091219.r'.
7. We have provided the output from the R CMD BATCH run of 'srac2lrac-example.R' in 'srac2lrac-example.Rout', along with the corresponding .RData file.
