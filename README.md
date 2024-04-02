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

For the RNA-AMP-Seq analysis, we adapted the IMPLICON pipeline avaiable here: https://github.com/FelixKrueger/IMPLICON

#### Overview of the pipeline

Assuming the input reads are
**R1:** `test_R1.fastq.gz` **R2:** `test_R2.fastq.gz`


**STEP I: UMI-Handling**

Read 2 has a 8nt UMI that has to be transfered to the read ID of both reads. Luckily, there is already a specialized command for this inside Trim Galore used in IMPLICON so we can just apply this:

```
trim_galore --paired --implicon *fastq.gz
```

**Output Files:**
```
test_8bp_UMI_R1.fastq.gz
test_8bp_UMI_R2.fastq.gz
```

**STEP II: Adapter-/quality trimming**

Standard Trim Galore usage to identify and remove read-through adapter contamination as well as poor quality base calls:

```
trim_galore --paired *UMI*fastq.gz
```

**Output Files:**
```
test_8bp_UMI_R1_val_1.fq.gz
test_8bp_UMI_R2_val_2.fq.gz
```

**STEP III: Genome Allignment**

Genome allignment was achieved with STAR instead of using Bismark as in the IMPLICON pipeline since Bismark is directed to then procede with the methylation exctraction, which is not relevant for this analysis.

```
STAR --runThreadN 8  \
        --genomeDir  /path/to/STAR/index \
        --outFileNamePrefix test_ \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesCommand zcat \
        --readFilesIn test_8bp_UMI_R1_val_1.fq.gz test_8bp_UMI_R2_val_2.fq.gz \
        --quantMode GeneCounts --sjdbGTFfile /path/to/annotation/file/gtf
```
**Output Files:**
```
test_Aligned.sortedByCoord.out.bam
```

**STEP IV: Deduplication using UMI-tools**

Since Bismark was not used for the allignments, we will use UMI-tools instead for the deduplication. 

First we need to index the allignemts. We can do so with SAM tools:
```
samtools index test_Aligned.sortedByCoord.out.bam
```

Then, we run UMI-tools:

```
umi_tools dedup --umi-separator=':' -I test_Aligned.sortedByCoord.out.bam --paired -S deduplicated_test_Aligned.sortedByCoord.out.bam
```

Finally, we need the .bam file to be indexed for the next step so we use sam tools again:
```
samtools index deduplicated_test_Aligned.sortedByCoord.out.bam
```

**Output Files:**
```
deduplicated_test_Aligned.sortedByCoord.out.bam
```

**STEP V: Allele Specific Counts**

[phASER](https://github.com/secastel/phaser/tree/master) was used as above described for the ASE analysis. In this example we use the ASD cell line (before and after cardiac differentiation): 

```
python phaser.py --vcf ASD.vcf.gz --bam deduplicated_ASD_Aligned.sortedByCoord.out.bam, \
deduplicated_ASDC_Aligned.sortedByCoord.out.bam --paired_end 1 --mapq 255 --baseq 10 --sample ASD --chr_prefix chr \
--blacklist hg38_hla.bed --haplo_count_blacklist hg38_haplo_count_blacklist.bed --threads 4 --o phaser_output_all/ASD
```

As before, we use the "phASER Gene AE" to produce gene-level haplotype expression quantifications.

```
python phaser_gene_ae.py --haplotypic_counts phaser_output_all/ASD.haplotypic_counts.txt \
--features human.gencode.v37.annotation.bed --o phaser_output_all/ASD.gene_ae.txt
```

For downstream analysis we discarded loci with a total read depth lower than 100 and the MAF was calculated as above stated.
