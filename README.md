# ErosionX
Repository dedicated to code and data associated with the manuscript "Erosion of X-Chromosome Inactivation in female hiPSCs is heterogeneous and persists during differentiation"

### Allele-Specific Expression (ASE) Analysis
Following the generation of a VCF file from WES data and the alignment of RNA-seq data to the reference human genome, we employed phASER for RNAseq-based phasing, to obtain gene-level haplotype expression data. We installed and ran phASER according to the usage tips here: https://github.com/secastel/phaser/tree/master/phaser

More specifically, to count allelic expression for isogenic control lines (F7cl15 and F7cl4 for instance), we use all bam files from those cell lines as input (--bam) and ‘F7’ as sample name (--sample) to isolate the variants of that cell line in the VCF file (--vcf).

```
python phaser.py --vcf hiPSCs.vcf.gz -–bam \
F7cl4_1.Aligned.sorted.bam, F7cl4_2.Aligned.sorted.bam, F7cl4_3.Aligned.sorted.bam, \
F7cl15_1.Aligned.sorted.bam,F7cl15_2.Aligned.sorted.bam,F7cl15_3.Aligned.sorted.bam \
--chr_prefix chr --hg38_hla.bed --haplo_count_blacklist blacklists/hg38_haplo_count_blacklist.bed \ 
--paired_end 1 --mapq 255 --baseq 10 --sample F7 --threads 8 --o phaser_out.F7
```

Authors suggest excluding variants in HLA genes using the "--blacklist" argument because of the high mapping error rate in these 
genes and exclude known problem sites from haplotypic counts using the "--haplo_count_blacklist" argument too. 
Location of these files can be found in their github page.

Next, we use the companion tool called "phASER Gene AE”, which takes the output files from phASER along with gene annotations
 and produces gene-level haplotype expression quantifications. Takes an input BED format file containing the coordinates for genes 
 where haplotypic counts are to be measured. Gene-level haplotypic counts are generated for each input BAM independently.

```
python phaser_gene_ae.py --haplotypic_counts phaser_out.F7.haplotypic_counts.txt \
–features human.gencode.v37.annotation.bed --o phaser_out.F7.gene_ae.txt
```

For downstream analysis we discarded loci with a total read depth lower than 10.
To quantify the effect size of allelic imbalance in expression for each gene, we computed the Minor Allele Frequency (MAF), 
which is computed as the ratio of minor allele read counts (the least common allele) to the total read counts from both alleles.
We defined genes with a MAF < 0.10 as fully monoallelic, genes with a MAF > 0.40 as fully biallelic, 
and the remaining genes as “in-between”. The final table is available as supplementary data.

### RNA AMPLICON-sequencing (RNA-AMP-seq) Analysis

For the RNA-AMP-Seq analysis, a simillar procedure was applied. A VCF file containing the SNPs of interest was used, which allowed to run the phaser command for all bam files using "All" as the sample argument (--sample). 

```
python phaser.py --vcf all_heterozygotic_VCF.vcf.gz --bam ASDC_S7_Aligned.sorted.bam, \
ASD_S1_Aligned.sorted.bam, CtrlDC_S10_Aligned.sorted.bam, CtrlD_S6_Aligned.sorted.bam, \
CtrlEC_S9_Aligned.sorted.bam, CtrlE_S5_Aligned.sorted.bam TCLABC_S8_Aligned.sorted.bam, \
TCLAB_S4_Aligned.sorted.bam --paired_end 1 --mapq 255 --baseq 10 --sample All --chr_prefix chr \
--blacklist hg38_hla.bed --haplo_count_blacklist hg38_haplo_count_blacklist.bed --threads 4 --o phaser_output_all/All
```

As before, we use the "phASER Gene AE" to produce gene-level haplotype expression quantifications.

```
python phaser_gene_ae.py --haplotypic_counts phaser_output_all/All.haplotypic_counts.txt \
--features human.gencode.v37.annotation.bed --o phaser_output_all/All.gene_ae.txt
```

For downstream analysis we discarded loci with a total read depth lower than 100 and the MAF was calculated as above stated.
