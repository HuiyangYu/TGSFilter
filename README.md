# HFKReads
<b> An ultra-fast and efficient tool for high frequency kmer sequencing reads extraction</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/HFKReads/releases/download/v2.02/HFKReads-2.02-Linux-x86_64.tar.gz
tar zvxf HFKReads-2.02-Linux-x86_64.tar.gz
cd HFKReads-2.02-Linux-x86_64
./hfkreads -h
```
### (2) Building from source （Linux or Mac）
```
git clone https://github.com/HuiyangYu/HFKReads.git
cd HFKReads
make
cd bin
./hfkreads -h
```
## 2. Usage
```
Usage: hfkreads -1 PE1.fq.gz -2 PE2.fq.gz -o OutFrefix
 Input/Output options:
   -1	<str>   paired-end fasta/q file1
   -2	<str>   paired-end fasta/q file2
   -s	<str>   single-end fasta/q
   -o	<str>   prefix of output file
 Filter options:
   -b	<int>   min base quality [0]
   -q	<int>   min average base quality [20]
   -l	<int>   min length of read [half]
   -r	<float> max unknown base (N) ratio [0.1]
   -k	<int>   kmer length [31]
   -w	<int>   window size [5]
   -m	<int>   min kmer count for high freq kmer [3] 
   -x	<int>   min count of read with high freq kmer [5]
   -n	<int>   read number to use [1000000]
   -a	        use all the read number
 Other options:
   -d  <int>    mode for deduplication (0:NO, 1:PE, 2:SE) [0]
   -c           compress the outPut File
   -f           outPut the KmerFre File
   -A           keep output quality info
   -t           thread to run [4]
   -h           show help [v2.02]
```
## 3. Example
### 3.1 Extracting high-frequency k-mer reads from paired-end sequencing reads
```
hfkreads -1 PE_1.fq.gz -2 PE_2.fq.gz -o test1
```
The output files consist of four files: <br>test1_pe_1.fa <br>test1_pe_2.fa <br>test1_se_1.fa <br>test1_se_2.fa <br>The output file with the 'se' label are unpaired high-frequency k-mer reads.

### 3.2 Extracting high-frequency k-mer reads from single-end sequencing reads
```
hfkreads -s SE.fq.gz -o test2
```
The output result consists of a single file named "test2.fa".

### 3.3 Extracting high-quality reads without filtering the k-mers
```
hfkreads -1 PE_1.fq.gz -2 PE_2.fq.gz -m 1 -o test3
hfkreads -s SE.fq.gz -m 1 -o test4
```
To filter out low-quality reads only, the default '-m 2' parameter should change to '-m 1'.

### 3.4 Filter ont ultra-long reads
```
hfkreads -s ont.fq.gz -m 1 -q 7 -l 20000 -a -c -A -o ont.filtered
```
The output result consists of a single file named "ont.filtered.fq.gz".
### 3.5 Filter Pacbio HIFI reads
```
hfkreads -s HIFI.fq.gz -m 1 -q 20 -l 10000 -a -c -A -o HIFI.filtered
```
The output result consists of a single file named "HIFI.filtered.fq.gz".

### 3.6 Other options
The '-A' parameter outputs quality values in fastq format, while the default is to output files in fasta format.<br>
The '-c' parameter compresses the output file, and by default, the output file is uncompressed.

## 4. License
-------

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
