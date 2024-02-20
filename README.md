# Influenza genomic surveillance in Galicia

The main goal of genomic surveillance of influenza is to report the geographical spread and temporal patterns of influenza clades and to detect resistance mutations. The objective is to contribute to reducing the burden of the seasonal influenza.

Surveillance of respiratory viruses has a multidisciplinary character, and must integrate clinical, epidemiological and virological data in an organized way to report to the ECDC by TESSY (The European Surveillance System).

As a part of the RELECOV (Red de Laboratorios de secuenciación de SARS-CoV-2) and for coordination with CNM-ISCIII (Centro Nacional de Microbiología, Instituto de Salud Carlos III), several activities are part of the routine of a sequencing laboratory: 
1. Sequencing and bioinformatic analysis
2. Submission sequences to GISAID
3. Report of clades to autonomous community Health Authorities (SERGAS)

Here we present an integrated influenza virus analysis pipeline including all these steps.

This is the result of a continuous collaborative effort of the following Institutions of the OMIC-G (Red de Laboratorios para la aplicación de Ómicas a la Microbiología Clínica en Galicia)

## Dependencies

The following software was used for this pipeline.

* BBSplit from [BBMap](https://sourceforge.net/projects/bbmap/)
* [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
* [samtools](https://www.htslib.org/)
* [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html)
* [nextclade](https://github.com/nextstrain/nextclade)
* [multiqc](https://multiqc.info/)
* [fluCLI](https://gisaid.org/)


## Sample files and folder structure

The starting point of this analysis are .fastq files stored in the folder *fastq*

The following folders are used through the whole analysis:

* *fastq*: with the original .fastq.gz files
* *gisaid*: for uploading the final .fasta and .csv files to GISAID
* *mid*: temporary files created in the process and .fa files of each segment
* *out*: final files created in the analysis: report, fasta, csv, ...
* *refs*: .fasta reference files for each segment


## 1. Flu typing: H1N1, H3N2 or BVic

We use the `bbsplit` utility from _BBmap_ in order to align the fastq sequences against each one of the HA segment of the 3 references (H1N1, H3N2 and BVic)

```
FILE="sample_name"
bbsplit.sh build=1 \
    in=fastq/${FILE}_L001_R1_001.fastq.gz \
    in2=fastq/${FILE}_L001_R2_001.fastq.gz \
    ref_1=refs/H1N1_HA.fasta \
    ref_2=refs/H3N2_HA.fasta \
    ref_3=refs/BVic_HA.fasta \
    out_1=fastq/pre/H1N1_HA_${FILE}.fastq.gz \
    out_2=fastq/pre/H3N2_HA_${FILE}.fastq.gz \
    out_3=fastq/pre/BVic_HA_${FILE}.fastq.gz
```

As a result, we get one .fastq.gz file for each sample and each type, stored in the folder _fastq/pre_ and named like this:
`${TYPE}_${SEGMENT}_${SAMPLE}.fastq.gz`

After that, we align the sequences in each obtained fast.gz file with `bwa-mem2` and study the coverage of each one with `samtools`. We only keep the ones with a coverage >80%. Usually one for each sample, except in case of coinfections.

```
FILE=${TYPE}_${SEGMENT}_${SAMPLE}

# Aligning:
bwa-mem2 mem -t 8 refs/${TYPE}_HA fastq/pre/${FILE}.fastq.gz > mid/${FILE}.sam

# 2. Converting from .sam to .bam:
samtools view -bS mid/${FILE}.sam -o mid/${FILE}.bam

# 3. Sorting:
samtools sort mid/${FILE}.bam -o mid/${FILE}.sorted.bam

# 4. Obtaining the coverage for H1N1/H3N2/BVic:
samtools index mid/${FILE}.sorted.bam
cover=$(samtools coverage -H mid/${FILE}.sorted.bam | cut -f6)

# 7. Discarding alignments below 80% coverage:
if (( $(echo "$cover < 80" |bc -l) )); then
    rm mid/${FILE}*
else
    rm mid/${FILE}.sam mid/${FILE}.bam
fi
```

As a result, in the _mid_ folder we get all the .bam files of the samples with this name structure:
`${TYPE}_HA_${SAMPLE}.sorted.bam`

## 2. Alignment against reference

Once we get the influenza subtype/s of each sample, we align the original fastq to the reference/s with coverage above 80%.

```
# In a loop through the mid folder, for each file "F" it's done the following:

# Split:
FILE=${F%.sorted.bam}
FILE=${FILE#mid/}
TYPE=${FILE:0:4}
SAMPLE=${FILE:8}

bbsplit.sh build=1 \
    in=fastq/${SAMPLE}_L001_R1_001.fastq.gz in2=fastq/${SAMPLE}_L001_R2_001.fastq.gz \
    ref_1=refs/${TYPE}_HA.fasta out_1=fastq/seg/${TYPE}_HA_${SAMPLE}.fastq.gz \
    ref_2=refs/${TYPE}_MP.fasta out_2=fastq/seg/${TYPE}_MP_${SAMPLE}.fastq.gz \
    ref_3=refs/${TYPE}_NA.fasta out_3=fastq/seg/${TYPE}_NA_${SAMPLE}.fastq.gz \
    ref_4=refs/${TYPE}_NP.fasta out_4=fastq/seg/${TYPE}_NP_${SAMPLE}.fastq.gz \
    ref_5=refs/${TYPE}_NS.fasta out_5=fastq/seg/${TYPE}_NS_${SAMPLE}.fastq.gz \
    ref_6=refs/${TYPE}_PA.fasta out_6=fastq/seg/${TYPE}_PA_${SAMPLE}.fastq.gz \
    ref_7=refs/${TYPE}_B1.fasta out_7=fastq/seg/${TYPE}_B1_${SAMPLE}.fastq.gz \
    ref_8=refs/${TYPE}_B2.fasta out_8=fastq/seg/${TYPE}_B2_${SAMPLE}.fastq.gz 

# From this point, each .fastq.gz file in the fastq/seg/ folder has the following name structure:
#   FILE={TYPE}_{SEGMENT}_{SAMPLE}

# Aligning:
bwa-mem2 mem -t 8 refs/${TYPE}_${SEGMENT} fastq/seg/${FILE}.fastq.gz > mid/${FILE}.sam

# Converting from .sam to .bam:
samtools view -bS mid/${FILE}.sam -o mid/${FILE}.bam

# Sorting:
samtools sort mid/${FILE}.bam -o mid/${FILE}.sorted.bam

# Triming primers:
ivar trim -i mid/${FILE}.sorted.bam -b refs/primers/${TYPE}_${SEGMENT}.primer.bed \
    -p mid/${FILE}.trimmed -e -m 32

# Re-sorting:
samtools sort mid/${FILE}.trimmed.bam -o mid/${FILE}.sorted.bam
```

As a result, we obtained, in the _mid_ folder one (or two, in coinfections) .sorted.bam file for each segment and each sample. Each file named like this:
`${TYPE}_${SEGMENT}_${SAMPLE}.sorted.bam`


## 3. Generation of .fasta consensus of each segment

We use `samtools` and `ivar` in order to get the .fasta file of the consensus of each segment of each sample, using the previously obtained .sorted.bam files

```
FILE=${TYPE}_${SEGMENT}_${SAMPLE}
samtools mpileup -aa -A -d 0 -Q 0 mid/${FILE}.sorted.bam \
    | ivar consensus -t 0 -p mid/${FILE} -i ${FILE}
```

As a result, in the _mid_ folder we get the .fa files with the same name structure:
`${TYPE}_${SEGMENT}_${SAMPLE}.fa`

Besides, we create a single .fasta file for each subtype with the HA segments of all samples. We'll need these files later to query nextclade. 

```
cat mid/H1N1_HA_*.fa > out/H1N1_HA_${FECHA}.fasta
cat mid/H3N2_HA_*.fa > out/H3N2_HA_${FECHA}.fasta
cat mid/BVic_HA_*.fa > out/BVic_HA_${FECHA}.fasta
```

## 4. Sumarize quality and coverage data

We obtain coverage data at 10X and 100X for each sample and store the results in a .csv file called _coverage.csv_. Later, we will use this file as part of the final report.

```
FILE=${TYPE}_${SEGMENT}_${SAMPLE}
samtools coverage -H mid/${FILE}.sorted.bam
samtools mpileup mid/${FILE}.sorted.bam
```


## 5. Geting information of each sample from nextclade

We use the `nextclade CLI` from Nextstrain to obtain data for each sample. For this, we need to provide the .fasta file of all the samples of each type, so we query Nextstrain for each different reference

```
nextclade run -d flu_h1n1pdm_ha -c out/nextclade_H1N1.csv out/H1N1_HA.fasta
nextclade run -d flu_h3n2_ha -c out/nextclade_H3N2.csv out/H3N2_HA.fasta
nextclade run -d flu_vic_ha -c out/nextclade_BVic.csv out/BVic_HA.fasta
```

The results, stored as .csv files in the _out_ folder, also will be used as a part of the final report.


## 6. Variant calling

We get data for minority variants, so we can make histogram plots in the final report based on this data.

```
samtools mpileup -A -d 0 --reference refs/${TYPE}_${SEGMENT}.fasta \
    -Q 0 mid/${FILE}.sorted.bam | ivar variants -p out/tsv/${FILE} \
    -t 0.03 -r refs/${TYPE}_${SEGMENT}.fasta -m 10 
```


## 7. Final report

The final report is made with an R script, summarizing all data obtained previously plus quality data obtained with `FastQC`, `Qualimap` and `multiqc` tools.

```
samtools stats ${FILE} > mid/$FECHA/${MUESTRA}.stats

fastqc fastq/*.fastq.gz --outdir=mid/fastqc/

qualimap multi-bamqc -d mid/qualimap_list.txt \
    -gff rsv_a/genemap.gff -outdir out -r mid/*.sorted.bam
    
multiqc --force -o out -n "multiqc.html" \
    -i "Report" -b "<a href='report.html'>Global Report</a>" mid 
```


## Flu sequences submission to GISAID

The objective of this part is to ease the process to send the metadata and .fasta segments to *GISAID* database. It's done with an R script and with the functionality of the `fluCLI` application from GISAID. 

In order to use the `fluCLI` command line utility, a GISAID user and a client_id it's needed.


### Metadata and fasta files

You need to complete the `gisaid_flu_template.csv` with the metadata of the samples (date, patient age, gender, Lab ID, ...). In our case, we use an .ods file from the LIS and process it with R to fill the template.csv.

The .fasta file with the 8 segments sequence of each sample is done by concatenating the individual .fa files in the _mid_ folder. The names of each sequence in the final .fasta file, should match the segment names in the template.csv. Again, we use R to do the work.


### Uploading to GISAID

The files generated in the previous step (.fasta and .csv) are uploaded to GISAID with the `fluCLI` utility. 

```
fluCLI upload --username XXXX --password YYYY --clientid ZZZZ \
    --metadata metadata.csv --fasta segments.fasta \
    --dateformat YYYYMMDD --log result.log
```

In the _result.log_ generated you'll have the accession_id of each sample and each segment.

After uploading to GISAID, you should access the GISAID website in order to release the samples.


## Final report to SAÚDE PÚBLICA DE GALICIA

In this final step, an .xlsx file is created with the data requested by SAÚDE PÚBLICA DE GALICIA. 
Again, we use R to write the .xlsx file including:
- accession_id numbers from the fluCLI log
- patient data and clade from the .ods from the LIS
- GISAID sample names from the .csv sent to GISAID


