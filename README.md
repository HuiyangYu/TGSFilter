# TGSFilter
<b> An ultra-fast and efficient tool for long reads filtering and trimming</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/TGSFilter/releases/download/v1.03/TGSFilter-1.013-Linux-x86_64.tar.gz
tar zvxf TGSFilter-1.03-Linux-x86_64.tar.gz
cd TGSFilter-1.03-Linux-x86_64
./hfkreads -h
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
Usage: tgsfilter -1 TGS_reads.fq.gz -o OutFile.fq.gz
 Options:
   -i	<str>   input of fasta/q file
   -o	<str>   output file
   -q	<int>   min Phred average quality score [10]
   -l	<int>   min length of read [100]
   -s	<int>   Trim N nucleotides from the start of a read [0]
   -e	<int>   Trim N nucleotides from the end of a read [0]
   -t		number of threads [1]
   -h		show help [v1.03]
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

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
