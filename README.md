# MetaRanker
A Computational Pipeline for Ranking Resistome Risk of Metagenomic Samples
=======

Workflow:
---------------
![workflow](https://github.com/user-attachments/assets/8d7a6b87-4f3f-4099-848c-ac5365bf72fe)

Requirements:
---------------
MetaRanker is written in python3. Please install these packages in python environment.
 - numpy
 - pandas
 - biopython

The following tools are needed for sequence processing. Please install them and add them to `$PATH`, so that these commands can be called in shell.
If you are using a Ubuntu operating system, these tools can be installed via `apt` (except for cd-hit).
 - [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
 - [cd-hit](https://github.com/weizhongli/cdhit)
 - [bwa](https://github.com/lh3/bwa)
 - [minimap2](https://github.com/lh3/minimap2) (for long reads)
 - [samtools](https://github.com/samtools/samtools)
 - pigz (recommend, for accelerating parsing gzipped fastq file)

Installation:
---------------
1.Clone this repository to your local directory. 

2.Install python packages. Skip if you have installed these packages.
```sh
pip install numpy pandas biopython
```
3.Install sequence processing tools.

4.Change directory to where you cloned this repository to.
```sh
cd /path/of/MetaRanker #example path
```
5.Uncompress `database.tar.gz` under current directory.
```sh
tar -zxvf database.tar.gz
```
6.Run python command to check installation and see help.
```sh
python MetaRanker.py -h
```

Usage: 
---------------
Basic usage:
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 32
```
For pair-ended reads:
```sh
python MetaRanker.py -c contigs.fa -r reads_1.clean.fastq.gz -R reads_2.clean.fastq.gz -o output_dir -t 32
```
For nanopore reads:
```sh
python MetaRanker.py -c contigs.fa -r reads.nanopore.clean.fastq.gz --nanopore -o output_dir -t 32
```
For pacbio reads:
```sh
python MetaRanker.py -c contigs.fa -r reads.pacbio.clean.fastq.gz --pacbio -o output_dir -t 32
```
Set minimum amount (1000) and cut-off length (250) of contigs:
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 32 --minnum 1000 --minlen 250
```
Keep original names of contigs. Usually names of contigs produced from assemblers contains white spaces,
which blastn will cut them off. So MetaRanker renames contigs by default.
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 32 --no_rename_contigs
```
Set weight of ARG, MGE, VF database to 1000, 1000, 1000:
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 32 --weight 1000 1000 1000
```

MetaRanker mainly use `subprocess` to call sequecne processing commands.
Some parameters used in MetaRanker can be passed, including `--blast_evalue`, `--blast_identity`,
`--blast_cover_len`, `--cdhit_identity` and `-t`/`--threads`

> [!NOTE]
> 1. Preprocess of reads
> To ensure the precision of sequence alignment, quality control of raw reads is needed.
> We recommend [fastp](https://github.com/OpenGene/fastp) or [trimmomatic](https://github.com/usadellab/Trimmomatic) for quality control of Illumina reads,
> and [chopper](https://github.com/wdecoster/chopper) for nanopore or pacbio reads.
> Also, host sequence removing can be applied if needed.
> 2. Assembly of contigs
> We recommend [megahit]() or [metaspades]() for assembling Illumina reads (megahit may be faster),
> and [flye]() for correcting, assembling and polishing 
Publications
---------------
This project was not published yet —— but you can still have a try on your metagenomic sequencing data.

License
---------------
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
