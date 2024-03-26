# TGSFilter
<b> An ultra-fast and efficient tool for long reads filtering and trimming</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/TGSFilter/releases/download/v1.08/TGSFilter-1.08-Linux-x86_64.tar.gz
tar zvxf TGSFilter-1.08-Linux-x86_64.tar.gz
cd TGSFilter-1.08-Linux-x86_64
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
Usage: tgsfilter -i TGS.raw.fq.gz -o TGS.clean.fq.gz
 Input/Output options:
   -i	<str>   input of fasta/q file
   -o	<str>   output of fasta/q file
 Basic filter options:
   -l	<int>   min length of read to out [1000]
   -L	<int>   max length of read to out
   -q	<int>   min mean base quality [auto]
   -n	<int>   read number for base content check [200000]
   -e	<int>   read end length for base content check [100]
   -5	<int>   trim bases from the front (5') of the read [auto]
   -3	<int>   trim bases from the tail (3') of the read [auto]
 Adapter filter options:
   -a	<str>   adapter sequence file 
   -A           disable reads filter, only for adapter identify
   -N	<int>   read number for adapter identify [200000]
   -E	<int>   read end length for adapter identify [100]
   -k	<int>   kmer size for adapter assembly [19]
   -y	<int>   min assembly adapter length [30]
   -m	<int>   min match length between read and adapter [4]
   -M	<int>   min match length between read middle and adapter [35]
   -s  <float>  min similarity between read end and adapter [0.75]
   -S  <float>  min similarity between read middle and adapter [0.9]
 Other options:
   -t           number of threads [3]
   -h           show help [v1.08]
```
## 3. Example

### 3.1 Filter HIFI reads
```
tgsfilter -i hifi.raw.fq.gz -o hifi.clean.fq.gz -t 10
```
The output files will be generated, namely 'hifi.clean.fq.gz'. <br>
The log information printed on the screen will contain the following content:
```
INFO: searching 5' adapter...
INFO: searching 3' adapter...
INFO: 5' adapter: ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT
INFO: 3' adapter: ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT
INFO: mean depth of 5' adapter: 18.4444
INFO: mean depth of 3' adapter: 14.4889
INFO: trim 5' end length: 150
INFO: trim 3' end length: 150
INFO: min mean base quality was set to: 20
INFO: min output reads length: 5000
10841385 reads with 199243618211 bases were input
164234 reads with 666540644 bases were dropped before filter
10677151 reads were trimmed 3203145300 bases in end
0 reads with 0 bases were dropped with low quality
3 reads were discard with 106304 bases with middle adapter
0 reads were trimmed 0 bases with adapter
0 reads with 0 bases were dropped before output
10677148 reads with 195373825963 bases were output
```
### 3.2 Filter ONT reads
```
tgsfilter -i ont.raw.fq.gz -o ont.clean.fq.gz -t 10
```
The log information printed on the screen will contain the following content:
```
INFO: searching 5' adapter...
INFO: searching 3' adapter...
INFO: 5' adapter: TGAAGCGGCGCACGAAAAACGCGAAAGCGTTTCACGATAAATGCGAAAAC
INFO: 3' adapter: 
INFO: mean depth of 5' adapter: 4686.72
INFO: mean depth of 3' adapter: 0
INFO: trim front length: 98
INFO: trim tail length: 12
INFO: min mean base quality was set to: 7
5301056 reads with 166540954904 bases were input
226266 reads with 199769142 bases were dropped before filter
5074790 reads were trimmed 558226900 bases in end
1053724 reads were trimmed 41307383 bases with adapter
358766 reads with 10869031124 bases were dropped with low quality
679963 reads with 17490611 bases were dropped before output
4707984 reads with 154855129744 bases were output
```
## 4. FAQ
### 4.1 How to set min mean base quality (-q)?
TGSFilter calculates the average quality value (meanQ) of the first 200,000 reads. If meanQ is greater than or equal to 25, then -q is set to 20; if meanQ is greater than or equal to 15 and less than 25, then -q is set to 10; otherwise, it is set to 7.<br>

The parameter '-q 20' is usually used to filter HiFi reads, while '-q 7' is typically used to filter ONT reads.<br>

### 4.2 How does TGSFilter identify adapter sequences?
TGSFilter employs three modes for identifying adapter sequences.<br>

(1) Users specify the adapter sequence through the '-a' parameter, which should be in fasta format.<br>

(2) Extract sequences of 100bp from the 5' and 3' ends, respectively, and align them with the general adapter library. If the global alignment similarity is over 90%, the adapter is detected.<br>

(3) Calculate the k-mer frequency of the extracted reads, filter out noisy k-mers, obtain candidate adapter sequences, and select the sequence with the highest k-mer depth as the adapter sequence.<br>

These three modes are executed sequentially. <br>

If you only want to identify adapters in your reads, you can add the '-A' parameter, which will instruct the program to only perform adapter search or assembly, taking less than 1 minute to complete.

### 4.3 How to set parameters for adapter assembly?
Due to the uncertainty of the quality and content of adapters in sequencing reads, although the quality of the adapter region is typically low, and the adapter content in HIFI reads is much lower than in ONT reads, excessively large or small k-mers are not suitable. In our tests, k-mers of 19 or 21 can successfully assemble the adapter sequences in both HIFI and ONT sequences. <br>

The read end length is also a factor affecting adapter identification, most adapters are present in the first 100bp of the 5' end and the last 100bp of the 3' end. Setting a detection length that is too large may mistakenly identify transposon sequences as adapters.
### 4.4 How to set the output file format? 
TGSFilter determines the output file format by recognizing the suffix of the output file name. For example:<br>

If the suffix of the output file name is '.fasta.gz' or '.fa.gz', a compressed fasta format file will be output. <br>

If the suffix of the output file name is '.fq' or '.fastq', an uncompressed fastq format file will be output.<br>
### 4.5 Why is the number of output reads greater than the number of input reads? 
In third-generation sequencing, adapters may be present in the middle of reads. When an adapter is detected in this region, we trim the read at that region, potentially resulting in an increased number of output reads compared to the input.

## 5. License
-------

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details
