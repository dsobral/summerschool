
# Summer School 2018

## High Throughput Sequencing (HTS) Data Analysis

### Practical prelude

**TASK**: Go to [http://cloud118.ncg.ingrid.pt](http://cloud118.ncg.ingrid.pt). To login use as *username* and *password* 'randomX', where X is the group number that was attributed to you (1,2,3....). 

For this session, we're going to use [Galaxy](https://galaxyproject.org/). This is a popular tool that allows even non-expert users to execute many different bioinformatics programs, particularly High Throughput Sequencing data analysis, through a simple web interface. Most software for the analysis of HTS data is freely available to users. Nonetheless, they often require the use of the command line in a Unix-like environment. User-friendly desktop software such as [CLC](https://www.qiagenbioinformatics.com/products/clc-genomics-workbench/) or [Ugene](http://ugene.net/) is available, but given the quick pace of developmpent in this area, they are constantly outdated. Moreover, even with better algorithms, HTS analysis must often be run in external servers due to intense computational requirements. Galaxy hides all the computational infrastructure from the end user, making it easier for researchers with limited training in informatics.

![Galaxy](galaxy.jpg)


**NOTE**: We're using a Galaxy instance specific for this session that is hosted in the [Portuguese Research Cloud (INCD)](http://www.incd.pt/) that will not be available after the course. Instances can be created according to research needs. There are also generic Galaxy instances publicly available such as [https://usegalaxy.eu/](https://usegalaxy.eu/). 


### Introduction

[Sanger sequencing](https://en.wikipedia.org/wiki/Sanger_sequencing) brought about a technological revolution, as it allowed to directly read DNA molecules with relative ease and affordability. The [Human Genome Project](https://en.wikipedia.org/wiki/Human_Genome_Project) motivated further progress, leading to automated DNA-sequencing machines capable of sequencing up to 384 samples in a single batch using capillary electrophoresis. 

At the moment, the high throughput sequencing technology most often used (by far) is Illumina. Similarly to the Sanger method, it is also based on the addition of nucleotides specifically modified to block DNA strand elongation, where each nucleotide is marked with a different color. Unlike the Sanger method, where a single DNA molecule is "read" at a time, modern illumina machines allow reading up to millions of DNA molecules simultaneously.  

![NGS Workflow](NGSworkflow.jpg)


The following links are a good source of information regarding this sequencing technology:
* [Illumina Sequencing by Synthesis](https://www.youtube.com/watch?&v=fCd6B5HRaZ8).
* [Elaine Mardis talk on NGS](https://www.youtube.com/watch?v=v1DbcJD4Ry0).

### The FastQ file format

Results of sanger sequencing are usually [fasta files](https://en.wikipedia.org/wiki/FASTA_format) (obtained from processing chromatograms). Most high-throughput sequencing machines output [fastq files](https://en.wikipedia.org/wiki/FASTQ_format), the “de facto” current standard in HTS. Like fasta, fastq files are simply text files, but where each block of information (a sequenced DNA fragment, or read) in this format is encoded as 4 lines:

	@read_identifier
	read_sequence
	+ separator line
	base_qualities
	
For example, here you have 8 lines of a fastq file, corresponding to 2 sequences:

	@HWI-M01876:76:000000000-AF16W:1:1101:10853:1000 1:N:0:CGTGACAGAT
	NTGTACTTCATCCGAAACTCGTGCTCATCTCTGCTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTGAT
	+
	#8ABCFGGGFCEDCFGGGGGGGFFCGEFGGGGGGFGGGGGGGGDEFGGGGGGGGGGGGGGGGGFFFEGGGGGGGGF
	@HWI-M01876:76:000000000-AF16W:1:1101:16471:1000 1:N:0:CGTGAACTTG
	NTTCCAGATATTCGATGCATGTGCCGCTCCTGTCGGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTGAT
	+
	#8BCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGG

Each base has a quality character associated with it, representing how confidently the machine identified (called) the base. The probability of error per base is given as a [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score), calculated from an integer value (Q) derived from the quality character associated to the base. The probability of error is given by the Phred score using P(Q)=10^(-Q/10). Useful reference values of Q include:
* Q=10 - 90% accuracy (0.1 error)
* Q=20 - 99% accuracy (0.01 error)
* Q=30 - 99.9% accuracy (0.001 error)
* Q=40 - 99.99% accuracy (0.0001 error)

Although there's theoretically no limit, Q usually goes up to around 40 in recent illumina machines.


**TASK**: Inspect the contents of mut_lib1_R1.fq. 

![View Data](view_data.jpg)

**Question**: What is the sequence of the first read in mut_lib1_R1.fq? And where is the information regarding the qualities of the bases?
<details><summary>Click Here to see the answer</summary><p>
The sequence of the first read is in the 2nd line:  "TTGCTGCTCTGGTTGTCGACAACGGCTCTGGCATGTGCAAGGCCGGATTTGCCGGAGACGATGCTCCCCGCGCCGTCTTCCCATCGATTGTGGGACGTC". The information on its base qualities is in the 4th line: "CCCFFFFFHHGHHIJIIJIIJJJGIHEHJJJHIIJGJJJJIJIJIJEHIHHHHHFCDBBDDDDDDDDDDDBB@BDBDDBDCDDDDB?BDBDDBBBBDDB".
</p></details>


### RNA-Seq

In RNA-Seq experiments, we're sequencing the transcripts that are expressed in chosen tissues/cells. The most commonly used illumina sequencing machines cannot directly read RNA (only DNA), therefore after purification, RNA needs to be converted to cDNA, after which it can be read like genomic DNA.

For this practical session, we're going to use data from [Guilgur et al, 2014](https://elifesciences.org/content/3/e02181). In this Drosophila melanogaster dataset, we have two conditions (WT and fandango mutant), with two replicates for each. To make it quick to run, we have extracted data for a very limited set of genes. You should be able to see 4 fastq files in your history (two WT replicates: WT_lib2_R1.fq, WT_lib1_R1.fq and two mutant replicates: mut_lib1_R1.fq, mut_lib2_R1.fq). The goal of this study is to indentify genes that are differentially expressed between the WT and the mutant condition.

The basic steps for this analysis:

  * Alignment of reads to the Drosophila genome

  * Generate table of counts of genes

  * Differential Analysis tests



