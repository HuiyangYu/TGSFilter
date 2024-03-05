# TGSFilter
<b> An ultra-fast and efficient tool for long reads filtering and trimming</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/TGSFilter/releases/download/v1.03/TGSFilter-1.013-Linux-x86_64.tar.gz
tar zvxf TGSFilter-1.03-Linux-x86_64.tar.gz
cd TGSFilter-1.03-Linux-x86_64
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
 Options:
   -i	<str>   input of fasta/q file
   -o	<str>   output of fasta/q file
   -l	<int>   min length of read to out [1000]
   -q	<int>   min mean base quality [auto]
   -5	<int>   drop bases from the front (5') of a read [0]
   -3	<int>   drop bases from the tail (3') of a read [0]
   -w	<int>   windows szie to cut off low quality region [50]
   -t           number of threads [3]
   -h           show help [v1.06]
```
## 3. Example

### 3.1 Filter ont ultra-long reads
```
tgsfilter -i ont.fq.gz -q 7 -l 20000 -o ont.filtered.fq.gz
```
### 3.2 Filter Pacbio HIFI reads
```
tgsfilter -i HIFI.fq.gz -q 20 -l 10000 -o HIFI.filtered.fq.gz
```

## 4. License
-------

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details
