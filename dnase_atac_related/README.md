# Scripts related to DNase- and ATAC-seq Analysis

--------------------------------------------------------------------------------

1. allele_specific_dnase.py

  **Decription:**

  Extract reads over a given set of SNPs with reference and variant allele specified.
  (DNase-, ATAC- or ChIP-seq). Perform basic filtering. Count the number of reads per allele and test for
  allele specificity using a two-sided binomial test.

  In the current version Indels are ignored. Also note that this script is fairly greedy and tests everything passing the inital filters. This may not lead to meaningful results for very low allele counts, where it is not clear if a low number of reads five sufficient support for calling a heterogeneous genotype.

  **Input:**

  SNP file (--snps) in bed-like or vcf format. Required columns are chr, pos (start, end for bed), reference and variant base (the column at which the ref and var allele are listed can be specified)

  BAM file(s) (--bam) input one or more bam files. Reads from all files will be piled-up so make sure they are comparable (e.g. from the same individual and a comparable data type)

  **Output:**

  Outputs a bed-like format given the SNP information a p-value and FDR corrected q-value for the respective specified reference and variant alleles. Alleles, allele counts and p-values results for all found alleles are reported in an additional column. (Format: Ref_base:Ref_count:Var_base:Var_count:Var_pvalue:Other_base:Other_count:Other_pvalue ...)

  **Requires:**

  Requires python3 (not tested for python2.X), pysam, scipy, statsmodels.

  **On CBRG Cluster** use module load python/3.3.2 and run as: *python3 allele_specific_dnase.py*
