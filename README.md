# Methods

Priya, S., Burns, M.B., Ward, T. et al. Identification of shared and disease-specific host gene–microbiome associations across human diseases using multi-omic integration. Nat Microbiol 7, 780–795 (2022). https://doi.org/10.1038/s41564-022-01121-z

## Gene expression data

"For host gene expression data for each
disease cohort, we used the ‘biomaRt’ R package (version 2.37.4) to only keep
data for protein-coding genes. We filtered out lowly expressed genes to retain
genes that are expressed in at least half of the samples in each disease cohort.
We performed variance stabilizing transformation using the R package ‘DESeq2’
(version 1.14.1) on the filtered gene expression read count data. We filtered out
genes with low variance, using 25% quantile of variance across samples in each
disease cohort as cut-off. Performing these steps for RNA-seq data for each disease
cohort separately resulted in a unique host gene expression matrix per disease for
downstream analysis, including 12,513 genes in the CRC dataset, 11,985 genes in
IBD dataset and 12,429 genes in IBS dataset."

## Microbiome data

"We performed the following steps for
microbiome data from each disease cohort separately. First, sequences that
were classified as either having originated from Archaea, chloroplasts, known
contaminants originating from laboratory reagents or kits, and soil- or
water-associated environmental contaminants were removed from the OTU
table as described earlier. Next, we summarized the OTU table at the species (if
present), genus, family, order, class and phylum taxonomic levels, and performed
prevalence and abundance-based filtering to retain taxa found at 0.001 relative
abundance in at least 10% of the samples..."

"To account for compositionality effects in microbiome datasets, we tested
two different approaches for performing centered log ratio (CLR) transformation
on our taxonomic data for each disease: (1) we concatenated the summarized
taxa matrices (count data) into a combined taxa matrix, and then applied CLR
transform on the combined matrix, (2) we CLR-transformed each taxon rank, and
then concatenated them into a combined matrix...Hence, we adopted the first approach for transforming our taxonomic
data."
