# **A detailed guide to assessing genome assembly based on long-read sequencing data using Inspector**

---

Scenario show cases for genome assembly based on long-read sequencing data using Inspector.

Inspector package can be found at [https://github.com/ChongLab/Inspector](https://github.com/ChongLab/Inspector).

## **Dependency**
Dependencies for Inspector:

- Python
    - Python module: pysam ≥ 0.22.1
    - Python module: statsmodels ≥ 0.14.2
- minimap2 ≥ 2.28
- SAMtools ≥ 1.20

Dependencies for Inspector error correction module:

- flye (tested with version 2.9.5)

Optional Dependencies (For Sequence Read Archive (SRA) Data):

If you need to download FASTQ reads from NCBI’s Sequence Read Archive (SRA):

- sra-tools (tested with 3.1.1)

## **Installation**

To create an environment with mamba or conda (recommended):

```bash
mamba create --name ins inspector
mamba activate ins
inspector.py -h
```

Git install after installing all the dependencies.

```bash
git clone https://github.com/ChongLab/Inspector.git
export PATH=$PWD/Inspector/:$PATH
```

A subset of human genome assembly is available as testing dataset to validate successful installation. The contig_test.fa includes two contigs (1.4Mbp and 10Kbp). The read_test.fastq.gz includes ~60X PacBio HiFi reads belonging to these two contigs. There are 3 structural errors and 281 small-scale errors present in the testing dataset.

```bash
cd Inspector/
./inspector.py -c testdata/contig_test.fa -r testdata/read_test.fastq.gz -o test_out/ --datatype hifi 
./inspector-correct.py -i test_out/ --datatype pacbio-hifi
```

## **General usage**

```bash
inspector.py [-h] -c contig.fa -r raw_reads.fa -o output_dict/
  required arguments:
  --contig,-c           FASTA/FASTQ file containing contig sequences to be evaluated
  --read,-r             A list of FASTA/FASTQ files containing long read sequences

  optional arguments:
  -h, --help            Show this help message and exit
  --version             Show program's version number and exit
  --datatype,-d         Input read type. (clr, hifi, nanopore, mixed) [clr]
  --ref                 OPTIONAL reference genome in .fa format
  --thread,-t           Number of threads. [8]
  --min_contig_length   Minimal length for a contig to be evaluated [10000]
  --min_contig_length_assemblyerror    Minimal contig length for assembly error detection. [1000000]
  --pvalue              Maximal p-value for small-scale error identification [0.01 for HiFi, 0.05 for others]
  --skip_read_mapping   Skip the step of mapping reads to contig
  --skip_structural_error       Skip the step of identifying large structural errors
  --skip_structural_error_detect       Skip the step of detecting large structural errors
  --skip_base_error     Skip the step of identifying small-scale errors
  --skip_base_error_detect      Skip the step of detecting small-scale errors from pileup

inspector-correct.py [-h] -i inspector_out/ --datatype pacbio-raw 
  required arguments:
  --inspector,-i        Inspector evaluation directory with original file names
  --datatype            Type of read used for Inspector evaluation. Required for structural error correction
  --outpath,-o          Output directory
  --flyetimeout         Maximal runtime for local assembly with Flye
  --thread,-t           Number of threads
  --skip_structural     Do not correct structural errors. Local assembly will not be performed
  --skip_baseerror      Do not correct small-scale errors
  

```

## **Use cases in this manuscript**

To demonstrate Inspector, we prepared four datasets to cover three unique aspects of long-read genome assembly applications. 

## **Scenario 1: Reference-free diploid and haplotype-resolved genome assemblies evaluation**


### **Stage 1: Download both assembly and raw genome reads files.**

Scenario 1 assembly files:

- hap1.fa.gz
- hap2.fa.gz

```bash
cd /path/to/assembled_contig/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/53FEE631-4264-4627-8FB6-09D7364F4D3B--ASM-COMP/HG002/assemblies/hifiasm_v0.19.5/hic/HG002.hap1.fa.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/53FEE631-4264-4627-8FB6-09D7364F4D3B--ASM-COMP/HG002/assemblies/hifiasm_v0.19.5/hic/HG002.hap2.fa.gz
```

Scenario 1 Raw reads files:

- m64011_190830_220126.Q20.fastq.gz
- m64011_190901_095311.Q20.fastq.gz
- m64012_190920_173625.Q20.fastq.gz
- m64012_190921_234837.Q20.fastq.gz

```bash
cd /path/to/rawfastq/
wget https://storage.googleapis.com/brain-genomics/awcarroll/t2t/fastq/q20/m64011_190830_220126.Q20.fastq.gz
wget https://storage.googleapis.com/brain-genomics/awcarroll/t2t/fastq/q20/m64011_190901_095311.Q20.fastq.gz
wget https://storage.googleapis.com/brain-genomics/awcarroll/t2t/fastq/q20/m64012_190920_173625.Q20.fastq.gz
wget https://storage.googleapis.com/brain-genomics/awcarroll/t2t/fastq/q20/m64012_190921_234837.Q20.fastq.gz
```

### **Stage 2: Genome assembly evaluation using Inspector for each haplotype**

Scenario 1 settings:

- Raw reads type used for evaluation: HiFi
- Reference: no
- Error correction: yes

```bash
inspector.py -c /path/to/assembled_contig/HG002.hap1.fa \
  -r /path/to/rawfastq/*.fastq.gz -o HG002_HiFi_hap1_noref/ -t8 --datatype hifi \
  1> HG002_HiFi_hap1_noref.log 2>&1

inspector-correct.py -t8 -i HG002_HiFi_hap1_noref/ --datatype pacbio-hifi \
  -o HG002_HiFi_hap1_noref

inspector.py -c /path/to/assembled_contig/HG002.hap2.fa \
  -r /path/to/rawfastq/*.fastq.gz -o  HG002_HiFi_hap2_noref/ -t8 --datatype hifi \
  1> HG002_HiFi_hap2_noref.log 2>&1

inspector-correct.py -t8 -i HG002_HiFi_hap2_noref/ --datatype pacbio-hifi \
  -o HG002_HiFi_hap2_noref
```

## **Scenario 2: Reference-free HiFi primary genome assembly evaluation with error correction**

### **Stage 1: Download both assembly and raw genome reads files.**

Scenario 2 assembly file:

- HiFi.hifiasm-0.12.pri.fa.gz

```bash
cd /path/to/assembled_contig/
wget https://zenodo.org/records/4393631/files/HG00733.HiFi.hifiasm-0.12.pri.fa.gz
```

Scenario 2 raw reads run IDs from Sequence Read Archive:

- ERR3822935
- ERR3861382
- ERR3861383
- ERR3861384
- ERR3861385
- ERR3861386
- ERR3861387

Download SRR_Acc_List.txt (select all 7 runs) form [https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERX3831682&o=acc_s%3Aa&s=ERR3822935](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERX3831682&o=acc_s%3Aa&s=ERR3822935) 

```bash
cd /path/to/rawfastq/
while read sample; do fastq-dump  $sample ; done <SRR_Acc_List.txt
```

### **Stage 2: Genome assembly evaluation using Inspector**

Scenario 2 settings:

- Raw reads type used for evaluation: HiFi
- Reference: no
- Error correction: yes

```bash
inspector.py -c /path/to/assembled_contig/HG00733.HiFi.hifiasm-0.12.pri.fa \
  -r /path/to/ rawfastq/*.fastq -o HG00733_HiFi_noref/ -t8 --datatype hifi \
  1>HG00733_HiFi_noref.log 2>&1

inspector-correct.py -t8 -i HG00733_HiFi_noref/ --datatype pacbio-hifi -o HG00733_HiFi_noref/
```

## **Scenario 3: Reference-based ONT genome assembly evaluation with error correction**

**Stage 1: Download assembly, raw genome reads files, and genome reference file**

Scenario 3 assembly filename:

- ONT.Flye-2.4.2_Helen.pri.fa.gz

```bash
cd /path/to/assembled_contig/
wget https://zenodo.org/records/4393631/files/HG00733.ONT.Flye-2.4.2_Helen.pri.fa.gz
```

Scenario 3 raw reads file names:

- HG00733_1.fastq.gz
- HG00733_2.fastq.gz
- HG00733_3.fastq.gz

```bash
cd /path/to/rawfastq/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/raw_data/nanopore/HG00733_1.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/raw_data/nanopore/HG00733_2.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/raw_data/nanopore/HG00733_3.fastq.gz
```

Scenario 3 reference FASTA file (T2T-CHM13 v2.0/hs1):

- File name: hs1.fa.gz

```bash
cd  /path/to/reference_genome/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
```

### **Stage 2: Genome assembly evaluation using Inspector**

Scenario 3 settings:

- Raw reads type used for evaluation: ONT
- Reference: yes
- Error correction: yes

```bash
inspector.py -c /path/to/assembled_contig/HG00733.ONT.Flye-2.4.2_Helen.pri.fa \   
  -r /path/to/ rawfastq/ *.fastq.gz --ref /path/to/reference_genome/hs1.fa \    
  -o HG00733_ONT_ref -t8 --datatype nanopore 1>HG00733_ONT_ref.log 2>&1
  
inspector-correct.py -t8 -i HG00733_ONT_ref/ --datatype nano-raw -o HG00733_ONT_ref/
```

## **Scenario 4: Reference-only ONT genome assembly evaluation with error correction**

### **Stage 1: Download assembly, raw genome reads files, and genome reference file**

Scenario 4 is a reference-only evaluation, where the assembled contigs and the reference genome are the same as those in Scenario dataset 3.

### **Stage 2: Genome assembly evaluation using Inspector**

Scenario 4 settings:

- Raw reads type used for evaluation: none
- Reference: yes
- Error correction: no

```bash
inspector.py -c /path/to/assembled_contig/HG00733.ONT.Flye-2.4.2_Helen.pri.fa \    
  -r emptyfile --ref /path/to/reference_genome/hs1.fa -o HG00733_ONT_ref_only \    
  -t8 1>HG00733_ONT_ref_only.log 2>&1
```

## **Output file descriptions**

The primary functionality of Inspector will generate four folders and ten files. Here we describe the most important outputs.

- Summary_statistics–- this file includes reports on contig continuity statistics
- small_scale_error.bed – this file lists all small-scale errors identified in the assembly
- small_scale_error_ref.bed – this file lists all small-scale errors identified in the assembly when a reference genome is provided
- structural_error.bed – this file documents all structural errors identified in the assembly
- structural_errors_ref.bed – this file lists structural errors identified in the assembly when a reference genome is provided
- contig_corrected.fa – this is the corrected assembly after Inspector’s error correction process

The “summary_statistics” file presents 26 metrics for the reference-free approach, the detailed explaination can be checked in ([Example output/Table 1. Summary statistics for reference free evaluation.csv](https://github.com/yuweisong-uab/Inspector_protocal/blob/main/Example%20output/Table%201.%20Summary%20statistics%20for%20reference%20free%20evaluation.csv)). When a reference genome is provided, an additional 13 statistics are reported ([Example output/Table 2. Summary statistics for reference-based evaluation.csv](https://github.com/yuweisong-uab/Inspector_protocal/blob/main/Example%20output/Table%202.%20Summary%20statistics%20for%20reference-based%20evaluation.csv))

## Questions

If you have any questions, please leave message in Issues here or contact Yuwei Song: [songyw@uab.edu](mailto:songyw@uab.edu).
