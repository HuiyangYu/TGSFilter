# TGSFilter
<b> An ultra-fast and efficient tool for long reads filtering and trimming</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/TGSFilter/releases/download/v1.10/TGSFilter-1.10-Linux-x86_64.tar.gz
tar zvxf TGSFilter-1.10-Linux-x86_64.tar.gz
cd TGSFilter-1.10-Linux-x86_64
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
 Input/Output options:
   -i	<str>   input of fasta/q file
   -o	<str>   output of fasta/q file
 Basic filter options:
   -l	<int>   min length of read to out [1000]
   -L	<int>   max length of read to out
   -q	<int>   min mean base quality [auto]
   -5	<int>   trim bases from the 5' end of the read [0]
   -3	<int>   trim bases from the 3' end of the read [0]
 Adapter filter options:
   -a	<str>   adapter sequence file 
   -A           disable reads filter, only for adapter identify
   -N	<int>   read number for adapter identify [200000]
   -E	<int>   read end length for adapter identify [100]
   -k	<int>   kmer size for adapter assembly [19]
   -y	<int>   min assembly adapter length [30]
   -m	<int>   min match length between read end and adapter [4]
   -M	<int>   min match length between read middle and adapter [35]
   -s  <float>  min similarity between read end and adapter [0.75]
   -S  <float>  min similarity between read middle and adapter [0.9]
   -D           split reads with middle adapter instead of discard
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
INFO: trim 5' end length: 0
INFO: trim 3' end length: 0
INFO: min mean base quality was set to: 20
INFO: min output reads length: 1000
INFO: 10841385 reads with a total of 199243618211 bases were input.
INFO: 83 reads were discarded with 52022 bases before filtering.
INFO: 0 reads were trimmed by 0 bases at the end.
INFO: 0 reads were discarded with 0 bases due to low quality.
INFO: 3 read was discarded with 107204 bases due to a middle adapter.
INFO: 209263 reads were trimmed by 6878822 bases with an adapter.
INFO: 0 reads were discarded with 0 bases before output.
INFO: 10841299 reads with a total of 199236580163 bases were output
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
INFO: mean depth of 5' adapter: 10721.6
INFO: mean depth of 3' adapter: 0
INFO: trim 5' end length: 0
INFO: trim 3' end length: 0
INFO: min mean base quality was set to: 10
INFO: min output reads length: 1000
INFO: 5301056 reads with a total of 166540954904 bases were input.
INFO: 134711 reads were discarded with 103200391 bases before filtering.
INFO: 0 reads were trimmed by 0 bases at the end.
INFO: 644139 reads were discarded with 19209834250 bases due to low quality.
INFO: 411 read was discarded with 11178458 bases due to a middle adapter.
INFO: 663024 reads were trimmed by 50778375 bases with an adapter.
INFO: 4636 reads were discarded with 4456634 bases before output.
INFO: 4517159 reads with a total of 147161506796 bases were output
```
## 4. FAQ
### 4.1 How to set min mean base quality (-q)?
TGSFilter calculates the average quality value (meanQ) of the first 10,000 reads with a length greater than or equal to 5000 bp. If meanQ is greater than or equal to 25, then -q is set to 20; otherwise, it is set to 10.<br>

The parameter '-q 20' is usually used to filter HiFi reads, while '-q 10' is typically used to filter ONT reads.<br>

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
