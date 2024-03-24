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
Usage: tgsfilter -1 TGS.raw.fq.gz -o TGS.clean.fq.gz
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
Two output files will be generated, namely 'hifi.raw.adapter.fa' and 'hifi.clean.fq.gz'. <br>
The log information printed on the screen will contain the following content:
```
INFO: searching 5' adapter...
INFO: searching 3' adapter...
INFO: 5' adapter: ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT
INFO: 3' adapter: ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT
INFO: mean depth of 5' adapter: 18.4444
INFO: mean depth of 3' adapter: 14.4889
INFO: trim front length: 12
INFO: trim tail length: 11
INFO: min mean base quality was set to: 20
10841385 reads with 199243618211 bases were input
90 reads with 59119 bases were dropped before filter
10841295 reads were trimmed 249349785 bases in end
990743 reads were trimmed 32623438 bases with adapter
0 reads with 0 bases were dropped with low quality
329014 reads with 2551884 bases were dropped before output
10841402 reads with 198959033985 bases were output
```
### 3.2 Filter ONT reads
```
tgsfilter -i ont.raw.fq.gz -o ont.clean.fq.gz -t 10
```
The log information printed on the screen will contain the following content:
```
INFO: searching 5' adapter...
INFO: searching 3' adapter...
INFO: not found adapter
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

## 4. License
-------

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details
