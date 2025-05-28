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
 - [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
 - [cd-hit](https://github.com/weizhongli/cdhit)
 - [bwa](https://github.com/lh3/bwa)
 - [minimap2](https://github.com/lh3/minimap2) (for long reads)
 - [samtools](https://github.com/samtools/samtools)
 - pigz (recommend, for accelerating parsing gzipped fastq file)

Installation:
---------------
1. Clone this repository to your local directory. 

2. Install python packages. Skip if you have installed these packages.
```sh
pip install numpy pandas biopython
```
3. Install blastn, cd-hit, bwa/minimap2, samtools and pigz. The installation instruction can be found at links above.
If you are using a Ubuntu operating system, these tools can be installed via `apt` (except for cd-hit):
```sh
sudo apt-get install ncbi-blast+
sudo apt-get install bwa
sudo apt-get install minimap2
sudo apt-get install samtools
sudo apt-get install pigz
```
If you downloaded executable binary tools or compiled tools from source, Please add them to `$PATH`:
```sh
export PATH=$PATH:/home/software/bwa # add executable binary to $PATH, edit in ~/.bashrc
```
```sh
cd /home/software/cd-hit
make # complie
make install # add complied tool to /usr/local/bin
```
Check if the tools are installed. If is installed, output should not be blank:
```sh
which blastn
which cd-hit
which bwa
which minimap2
which samtools
which pigz
```
4. Change directory to where you cloned this repository to.
```sh
cd /path/of/MetaRanker #example path
```
5. Run python command to check installation and see help.
```sh
python MetaRanker.py -h
```

Usage: 
---------------
1. Basic usage:
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 16
```
2. For pair-ended reads:
```sh
python MetaRanker.py -c contigs.fa -r reads_1.clean.fastq.gz -R reads_2.clean.fastq.gz -o output_dir -t 16
```
3. For nanopore reads:
```sh
python MetaRanker.py -c contigs.fa -r reads.nanopore.clean.fastq.gz --nanopore -o output_dir -t 16
```
4. For pacbio reads:
```sh
python MetaRanker.py -c contigs.fa -r reads.pacbio.clean.fastq.gz --pacbio -o output_dir -t 16
```
5. Set minimum amount (1000) and cut-off length (250) of contigs:
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 16 --minnum 1000 --minlen 250
```
6. Keep original names of contigs. Usually names of contigs produced from assemblers contains white spaces,
which blastn will cut them off. So MetaRanker renames contigs by default.
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 16 --no_rename_contigs
```
7. Set weight of ARG, MGE, VF database to 1000, 1000, 1000:
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 16 --weight 1000 1000 1000
```
8. Pass parameters to sequence processing tools.
MetaRanker mainly use `subprocess` to call sequence processing commands.
Some parameters used in MetaRanker can be passed, including `--blast_evalue`, `--blast_identity`,
`--blast_cover_len`, `--cdhit_identity` and `-t`/`--threads`
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 64 --blast_evalue 1e-5 --blast_identity 0.9 --blast_cover_len 85 --cdhit_identity 0.9
```
9. Force to overwrite existing output files of blastn, cd-hit, bwa, minimap2 and samtools. 
By default, MetaRanker will not call a sequence processing command if the output file exists.
```sh
python MetaRanker.py -c contigs.fa -r reads.clean.fastq.gz -o output_dir -t 16 --force
```

> [!NOTE]
> 1. Preprocess of reads
> To ensure the precision of sequence alignment, quality control of raw reads is needed.
> We recommend [fastp](https://github.com/OpenGene/fastp) or [trimmomatic](https://github.com/usadellab/Trimmomatic) for quality control of Illumina reads,
> and [chopper](https://github.com/wdecoster/chopper) for nanopore or pacbio reads.
> Also, host sequence removing can be applied if needed.
> 2. Assembly of contigs
> We recommend [megahit](https://github.com/voutcn/megahit) or [metaspades](https://github.com/ablab/spades) for assembling Illumina reads (megahit may be faster),
> and [flye](https://github.com/mikolmogorov/Flye) for correcting, assembling and polishing nanopore or pacbio reads.

An example pipeline:
```sh
fastp -i sample_1.fastq.gz -o sample_1.clean.fastq.gz -I sample_2.fastq.gz -O sample_2.clean.fastq.gz -w 16
megahit -1 sample_1.clean.fastq.gz -2 sample_2.clean.fastq.gz -o contig_dir --out-prefix sample -t 16
python MetaRanker.py -c contig_dir/sample.contigs.fa -r sample_1.clean.fastq.gz -R sample_2.clean.fastq.gz -o metaranker_output -t 16
```

Outputs: 
---------------
In output directory, result files should be:

1. RiskVector, RiskModule, CoocurScore, RiskIndex, ReadsNum, BasesNum, ContigsNum of sample
> metaranker_output/risk_result/RiskStat_sample.tsv

2. Co-occurrence matrix of sample
> metaranker_output/risk_result/RiskMatrix_sample.csv

3. BLAST results in M8 format
> metaranker_output/output_M8/*.tsv

4. Filtered, annotated M8 tables
> metaranker_output/preprocessed_M8/*.tsv

5. Filtered, annotated M8 tables with categorized genes
> metaranker_output/preprocessed_M8/categorized_M8/*.tsv


6. Sequences and depths of Risk Elements
> metaranker_output/risk_elements/*

7. BPM abundance of Risk Elements
> metaranker_output/BPM/*.tsv

8. High risk sequences with co-ocurrence structures dumped from contigs
> metaranker_output/coocur_structures/*.fasta

9. Contigs renamed by MetaRanker (can be removed as needed)
> metaranker_output/temp/*.fa


> [!TIP]
> The co-ocurrence structures of each samples can be visualized with `SeqVisualize.py`.
> 
> After concatenating `RiskStat_*.tsv` of samples, a 3D hazard space plot can be produced with `Plot3Dspace.py`.

Publications
---------------
This project was not published yet —— but you can still have a try on your metagenomic sequencing data.

License
---------------
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
