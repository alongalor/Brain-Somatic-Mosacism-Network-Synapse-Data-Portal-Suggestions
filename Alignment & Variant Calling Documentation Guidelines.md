# Alignment and Variant Calling Documentation Guidelines:

Documentation should contain a minimal amount of information that allows others to quickly understand and easily reproduce pipelines. This should include all commands applied and software versions used as well as links to code for custom scripts.

Please see the Walsh-Park piplines below for examples:



## Alignment Pipeline

### Split Fastqs

Decompress and split fastq files `R1.fastq.gz` and `R2.fastq.gz` into multiple fastq files: `1.R1.fq, 2.R1.fq, ... , n.R1.fq` and `1.R2.fq, 2.R2.fq, ... , n.R2.fq`, respectively.

```
zcat R1.fastq.gz | python ali.py -step1_split_fq -readRno R1 -conf conf.txt

zcat R2.fastq.gz | python ali.py -step1_split_fq -readRno R2 -conf conf.txt
```

### For each pair of split fq files, `i.R1.fq, i.R2.fq`, for all `i` in `[1, 2, ... , n]`, perform the following:

##### BWA

```
bwa mem -M -t 2 ref_genome i.R1.fq i.R2.fq | samtools view -bSho i_bwa.bam -
```

##### Sort BAM 

```
java -Xmx4G -Xms4G -Djava.io.tmpdir=/tmp -jar picard SortSam SORT_ORDER=coordinate INPUT=i_bwa.bam OUTPUT=i_sort.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
```

##### Add Read Groups

```
java -jar -Xmx4G -Xms4G picard AddOrReplaceReadGroups RGID=L7 RGLB=ID RGPL=illumina RGPU=unit1 RGSM=TEST I=i_sort.bam O=i_addrg.bam

samtools index i_addrg.bam
```

### Mark Duplicates

```
java -Xmx20G -Xms20G -jar picard MarkDuplicates \
M=out.dup.matrix \
COMPRESSION_LEVEL=1 \
ASSUME_SORTED=True \
VALIDATION_STRINGENCY=LENIENT \
REMOVE_DUPLICATES=false \
I=1_addrg.bam \
I=2_addrg.bam \
...
I=n_addrg.bam \
O=mkdup.bam \
TMP_DIR=/tmp

samtools index mkdup.bam
```

Note that this step merges `i_addrg.bam` for all `i` is in `[1, 2, ... , n]` into `mkdup.bam`

### Realign Indels

```
java -Xmx4G -Xms4G -jar gatk \
-T RealignerTargetCreator \
-R ref_genome \
-I mkdup.bam \
-o indel.intervals \
 
java -Xmx4G -Xms4G -jar gatk \
-T IndelRealigner \
-R ref_genome \
-I out.bam \
-targetIntervals indel.intervals -o indel.bam

samtools index indel.bam
```

### Recalibrate Bases
```
java -Xmx20G -Xms20G -jar gatk \
-T BaseRecalibrator \
-R ref_genome \
-knownSites dbsnp  \
-I indel.bam \
-o bqsr.bqsr \
-nct 8

java -Xmx10G -Xms10G -jar gatk \
-T PrintReads \
-R ref_genome \
-BQSR bqsr.bqsr \
-I indel.bam \
-o bqsr.bam \
-nct 8
```

#### Software Versions:

*java*: jdk1.7.0_71

*bwa*: 0.7.8

*picard*: 1.13

*samtools*: 1.3

*gatk*: 3.5

*dbsnp*: 147

*ref_genome*: human_g1k_v37_decoy

#### Custom Scripts:

Please see `ali.py`, attached, for details on the "Split Fastqs" step.



## Varaint Calling Pipelines

### Tumor/Normal Variant Calling

```
java -Xms4G -Xmx4G -jar gatk \
-T MuTect2 \
-R ref_genome \
-I:normal normal.bam \
-I:tumor tumor.bam \
--cosmic cosmic_coding_mutations \
--dbsnp dbsnp \
-o tumor_normal.vcf \
-nct 2
```

### Tumor/Normal Variant Calling with Panel of Normals Creation

#### For each normal sample to be added to the Panel of Normals, `normal_i.bam` for all `i` in `[1, 2, ... , n]`, perform the following:

##### Call Normal Samples Separately

```
java -Xms4G -Xmx4G -jar gatk \
-T MuTect2 \
-R ref_genome \
-I:tumor normal_i.bam \
--dbsnp dbsnp \
--cosmic cosmic_coding_mutations \
--artifact_detection_mode \
-o pon.normal_i.vcf \
-nct 2 \
-dt NONE 
```

#### Combine Variants

```
java -jar GenomeAnalysisTK.jar \
 -T CombineVariants \
 -R reference.fasta \
 -V pon.normal_1.vcf -V pon.normal_2.vcf ... -V pon.normal_n.vcf \
 -minN 2 \
 --setKey "null" \
 --filteredAreUncalled \
 --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
 [-L targets.interval_list] \
 -o MuTect2_PON.vcf
```

Note that this step merges `pon.normal_i.vcf` for all `i` in `[1, 2, ... , n]` into `MuTect2_PON.vcf`

#### Call Tumor vs. Panel of Normals  

```
java -Xms4G -Xmx4G -jar gatk \
-T MuTect2 \
-R ref_genome \
--dbsnp dbsnp \
--cosmic cosmic_coding_mutations \
-I:tumor tumor.bam \
-PON MuTect2_PON.vcf \
-nct 2 \
-o tumor_normal.pon.vcf \
-dt NONE
```

### Software Versions

*java*: jdk1.7.0_71

*gatk*: 3.5

*cosmic_coding_mutations*: v72

*dbsnp*: v147