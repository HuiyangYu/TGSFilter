# TGSFilter
<b> An ultra-fast and efficient tool for long reads filtering and trimming</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/TGSFilter/releases/download/v1.11/TGSFilter-1.11-Linux-x86_64.tar.gz
tar zvxf TGSFilter-1.11-Linux-x86_64.tar.gz
cd TGSFilter-1.11-Linux-x86_64
./tgsfilter -h
```
### (2) Building from source （Linux or Mac）
```
git clone https://github.com/HuiyangYu/TGSFilter.git
cd TGSFilter
make
cd bin
./tgsfilter -h
```
## 2. Usage
```
Usage: tgsfilter -i TGS.raw.fq.gz -x ont -o TGS.clean.fq.gz
 Input/Output options:
   -i   <str>   input of bam/fasta/fastq file
   -x   <str>   read type (ont|clr|hifi)
   -o   <str>   output of fasta/fastq file instead of stdout
 Basic filter options:
   -l   <int>   min length of read to out [1000]
   -L   <int>   max length of read to out
   -q  <float>  min Phred average quality score
   -Q  <float>  max Phred average quality score
   -n   <int>   read number for base content check [100000]
   -e   <int>   read end length for base content check [150]
   -b  <float>  bias (%) of adjacent base content at read end [1]
   -5   <int>   trim bases from the 5' end of the read
   -3   <int>   trim bases from the 3' end of the read
 Adapter filter options:
   -a   <str>   adapter sequence file 
   -A           disable reads filter, only for adapter identify
   -N   <int>   read number for adapter identify [100000]
   -E   <int>   read end length for adapter trim [150]
   -m   <int>   min match length for end adapter [15]
   -M   <int>   min match length for middle adapter [35]
   -T   <int>   extra trim length for middle adpter on both side [50]
   -s  <float>  min similarity for end adapter
   -S  <float>  min similarity for middle adapter
   -D           discard reads with middle adapter instead of split
 Downsampling options:
   -g   <str>   genome size (k/m/g)
   -d   <int>   downsample to the desired coverage (requires -g) 
   -r   <int>   downsample to the desired number of reads 
   -R  <float>  downsample to the desired fraction of reads 
   -k   <int>   kmer size for repeat evaluations [11] 
   -p   <int>   min repeat length of reads [0] 
   -F           disable reads filter, only for downsampling
 Other options:
   --qc         disable all filter, only for quality control 
   -f           force FASTA output (discard quality) 
   -c   <int>   compression level (0-9) for compressed output [6]
   -t   <int>   number of threads [16]
   -h           show help [v1.11]
```
## 3. Example

### 3.1 Filter HIFI reads
```
tgsfilter -i hifi.fq.gz -o hifi.clean.fq.gz -t 16 -x hifi
```
Two output files will be generated, namely <b>'hifi.clean.fq.gz'</b> and <b>'[hifi.clean.html](http://htmlpreview.github.io/?https://github.com/HuiyangYu/TGSFilter/blob/main/html/hifi.clean.html)'</b>. <br>

The log information printed on the screen will contain the following content:
```
INFO: read type: PacBio highly accurate long reads (hifi).
INFO: trim 5' end length: 7
INFO: trim 3' end length: 8
INFO: min output reads length: 1000
INFO: base quality scoring: Phred33
INFO: min Phred average quality score: 20
INFO: 5' adapter: 
INFO: 3' adapter: 
INFO: mean depth of 5' adapter: 0
INFO: mean depth of 3' adapter: 0
INFO: set PacBio blunt adapter to trim: ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
INFO: 19053754 reads with a total of 353996424711 bases were input.
INFO: 0 reads were discarded with 0 bases due to low quality.
INFO: 3 reads have adapter at 5', 3' and middle.
INFO: 4 reads have adapter at 5' and middle.
INFO: 11 reads have adapter at 3' and middle.
INFO: 138 reads have adapter at 5' and 3' end.
INFO: 360 reads only have adapter at middle.
INFO: 51806 reads only have adapter at 5' end.
INFO: 49952 reads only have adapter at 3' end.
INFO: 18951480 reads didn't have any adapter.
INFO: 298375518 bases were trimmed due to the adapter or base content bias.
INFO: 189 reads were discarded with 111248 bases due to the short length before output.
INFO: 19053949 reads with a total of 353697937945 bases after filtering.
INFO: Filtered reads were written to: hifi.clean.fq.gz.
INFO: Quality control report was written to: hifi.clean.html
```
### 3.2 Filter ONT reads
```
tgsfilter -i ont.fq.gz -l 50000 -o ont.clean.fq.gz -t 16 -x ont
```
<b>Attention:</b> If the sequencing depth of your ONT reads is < 40 X, setting the minimum length to 50 kb might be too stringent. You can set ```-l``` to 20,000 or 30,000 instead.<br>

Two output files will be generated, namely <b>'ont.clean.fq.gz'</b> and <b>'[ont.clean.html](http://htmlpreview.github.io/?https://github.com/HuiyangYu/TGSFilter/blob/main/html/ont.clean.html)'</b>.

The log information printed on the screen will contain the following content:
```
INFO: read type: NanoPore reads (ont).
INFO: trim 5' end length: 96
INFO: trim 3' end length: 10
INFO: min output reads length: 50000
INFO: base quality scoring: Phred33
INFO: min Phred average quality score: 10
INFO: 5' adapter: GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA
INFO: 3' adapter: 
INFO: mean depth of 5' adapter: 55085.2
INFO: mean depth of 3' adapter: 0
INFO: 3953401 reads with a total of 177832219038 bases were input.
INFO: 188645 reads were discarded with 6626347642 bases due to low quality.
INFO: 1 reads have adapter at 5', 3' and middle.
INFO: 968 reads have adapter at 5' and middle.
INFO: 1 reads have adapter at 3' and middle.
INFO: 3160 reads have adapter at 5' and 3' end.
INFO: 262 reads only have adapter at middle.
INFO: 3211881 reads only have adapter at 5' end.
INFO: 580 reads only have adapter at 3' end.
INFO: 547903 reads didn't have any adapter.
INFO: 400695818 bases were trimmed due to the adapter or base content bias.
INFO: 2672989 reads were discarded with 48381650389 bases due to the short length before output.
INFO: 1092920 reads with a total of 122423525189 bases after filtering.
INFO: Filtered reads were written to: ont.clean.fq.gz.
INFO: Quality control report was written to: ont.clean.html
```
## 4. FAQ
### 4.1 How to set min reads quality (-q)?
TGSFilter allows you to set the minimum read quality using the '-q' option. By default, the '-q' option is set to 0. However, if you do not specify the '-q' parameter but specify the type of reads using the '-x' parameter, the software will automatically set the  minimum read quality:<br> <b>
'-q 10' for ONT reads,<br> 
'-q 20' for HiFi reads,<br> 
'-q 0' for CLR reads.</b>

### 4.2 How to set min reads length (-l)?
The default minimum read length outputted by TGSFilter is 1000 bp. Generally, HIFI reads do not require additional length filtering. However, CLR and ONT reads require length filtering based on the sequencing depth. For CLR reads, you can set '-l 10000'; for ONT, '-l' can be set between 20000 and 50000. If the sequencing depth for ONT ultra-long reads exceeds 100X, you can also set '-l 100000'.

### 4.3 Should I keep as many reads as possible?
In general, whether for variant calling or genome assembly, 40-50X of sequencing depth is usually sufficient. Excessive depth may adversely affect downstream analyses. Filtering out low-quality and short sequencing reads while ensuring adequate sequencing depth not only reduces the time required for downstream analysis but also yields higher quality analytical results.

### 4.2 How does TGSFilter identify adapter sequences?
TGSFilter employs three steps for identifying adapter sequences.<br>

(1) Check whether users specify the adapter sequence through the '-a' parameter, which should be in fasta format. <br>

(2) Extract sequences of 100bp from the 5' and 3' ends, respectively, and align them with the general adapter library. If the global alignment similarity is over 90%, the adapter is detected.<br>

(3) Calculate the k-mer frequency of the extracted reads, filter out noisy k-mers, obtain candidate adapter sequences, and select the sequence with the highest k-mer depth as the adapter sequence.<br>

These three steps are executed sequentially. If users specify the adapter sequence file through the '-a' parameter, the following two steps will not be executed; otherwise, only the following two steps will be performed. <br>

If you only want to identify adapters in your reads, you can add the '-A' parameter, which will instruct the program to only perform adapter search or assembly, taking less than 1 minute to complete.

### 4.3 How to set parameters for adapter assembly?
Due to the uncertainty of the quality and content of adapters in sequencing reads, although the quality of the adapter region is typically low, and the adapter content in HIFI reads is much lower than in ONT reads, excessively large or small k-mers are not suitable. In our tests, k-mers of 19 or 21 can successfully assemble the adapter sequences in both HIFI and ONT sequences. <br>

The read end length is also a factor affecting adapter identification, most adapters are present in the first 100bp of the 5' end and the last 100bp of the 3' end. Setting a detection length that is too large may mistakenly identify transposon sequences as adapters.
### 4.4 How to set the output file format? 
TGSFilter determines the output file format by recognizing the suffix of the output file name. For example:<br>

If the suffix of the output file name is '.fasta.gz' or '.fa.gz', a compressed fasta format file will be output. <br>

If the suffix of the output file name is '.fq' or '.fastq', an uncompressed fastq format file will be output.<br>

## 5. License
-------

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details
