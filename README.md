## NeuralPolish ##
NeuralPolish: a novel Nanopore polishing method based on alignment matrix construction and orthogonal Bi-GRU Networks.

## Overview ##
Here, we developed a novel polishing method, named NeuralPolish, based on alignment matrix construction and orthogonal Bi-GRU networks. In this method, we designed an alignment feature matrix for representing read-to-assembly alignment. Each row of the matrix represents a read, and each column represents the aligned bases at each position of the contig. In the network architecture, a bi-directional GRU network is used to extract the sequence information inside each read by processing the alignment matrix row by row. After that, the feature matrix is processed by another bi-directional GRU network column by column to calculate the probability distribution. Finally, a CTC decoder generates a polished sequence with a greedy algorithm.

## Installation ##
Using this method requires the user to install several tools:
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [Racon](https://github.com/isovic/racon)
- [GCC>=7.3](http://gcc.gnu.org/releases.html)

dependencies:
```
pip install pyyaml pysam python-Levenshtein numpy biopython tensorboardX
```

if the machine has GPUs, you can install pytorch-gpu >=1.4.0 environment:
```
conda install pytorch=1.4.0
```

if the machine only has Cpus, install pytorch-cpu >=1.4.0 environment:
```
conda install pytorch-cpu=1.4.0
```

Install BlockPolish from the GitHub repository:
```
git clone https://github.com/huangnengCSU/NeuralPolish.git
cd neuralpolishextract/ && cmake . && make && cd ..
```

## Useage ##
```
bash polish.sh -h
usage: ./polish.sh -i <draft_assembly> -r <basecalled_reads> -o <output_dir> [options ...]
options:
  -i STR   input file in FASTA format, containing draft assembly which will be polished
  -r STR   input file in FASTA/FASTQ format, containing basecalled raw reads used for polishing
  -o STR   output directory of files
  -q INT   skip alignments with mapQ smaller than INT, default: 0
  -Q INT   skip bases with baseQ/BAQ smaller than INT, default: 13
  -t INT   number of threads, default: 16
  -g INT   index of gpu device, default: None
  -c INT   coverage of reads used for polishing, default: 40
  -f STR   neural network config file, default: ./brnn/config/wtdbg2+racon.3x3.yaml
  -m STR   neural network pre-trained file, default: ./brnn/ctc.ecoli_yeast_chr21.wtdbg2+racon.layer3x3.model/model.chkpt
```

## Examples ##
Here, we use NeuralPolish to polish the draft assembly of *E. coli* assembled by Flye. The three essential parameters are draft assembly file (-i), basecalled raw reads (-r) and output file directory (-o). The parameter of gpu device index (-g) specifies which GPU device will be used to perform the computation of neural network prediction.
```
bash polish.sh -i ~/data/ecoli/ecoli_flye/assembly.fasta -r ~/data/ecoli/call.fastq -o polished_ecoli -g 0
```

## License ##
Copyright (C) 2020 by Huangneng (huangn94@foxmail.com)


