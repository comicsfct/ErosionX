# ErosionX
Repository dedicated to code and data associated with the manuscript "Erosion of X-Chromosome Inactivation in female hiPSCs is heterogeneous and persists during differentiation"

### Allele-Specific Expression (ASE) Analysis
Following the generation of a VCF file from WES data and the alignment of RNA-seq data to the reference human genome, we employed phASER for RNAseq-based phasing, to obtain gene-level haplotype expression data. We installed and ran phASER according to the usage tips here: https://github.com/secastel/phaser/tree/master/phaser

More specifically, to count allelic expression for isogenic control lines (F7cl15 and F7cl4 for instance), we use all bam files from those cell lines as input (--bam) and ‘F7’ as sample name (--sample) to isolate the variants of that cell line in the VCF file (--vcf).

```
python phaser.py --vcf hiPSCs.vcf.gz \
–bam F7cl4_1.Aligned.sortedByCoord.out.bam, F7cl4_2.Aligned.sortedByCoord.out.bam, F7cl4_3.Aligned.sortedByCoord.out.bam, \
F7cl15_1.Aligned.sortedByCoord.out.bam,F7cl15_2.Aligned.sortedByCoord.out.bam,F7cl15_3.Aligned.sortedByCoord.out.bam \
--chr_prefix chr --hg38_hla.bed --haplo_count_blacklist blacklists/hg38_haplo_count_blacklist.bed --paired_end 1 --mapq 255 --baseq 10 \
--sample F7 --threads 8 --o phaser_out.F7
```





