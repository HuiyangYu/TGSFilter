#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <ctime>
#include <thread>
#include <algorithm> //for std::sort
#include <cstdlib>  // for getenv
#include <unistd.h> // for access
#include <cstring> // for strlen and mmemcpy
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cstdio>
#include <map>
#include <random>
#include "sam.h"
#include "hts.h"
#include "zlib.h"
#include "edlib.cpp"
#include "igzip_lib.h"
#include "libdeflate.h"
#include "concurrentqueue.h"

using namespace std;

uint8_t Base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};

int  TGSFilter_usage() {
	cout <<""
		"Usage: tgsfilter -1 TGS.raw.fq.gz -o TGS.clean.fq.gz\n"
		" Input/Output options:\n"
		"   -i  <str>   input of bam/fasta/fastq file\n"
		"   -o  <str>   output of fasta/fastq file instead of stdout\n"
		" Basic filter options:\n"
		"   -l  <int>   min length of read to out [1000]\n"
		"   -L  <int>   max length of read to out\n"
		"   -q  <int>   min Phred average quality score [0]\n"
		"   -Q  <int>   max Phred average quality score \n"
		"   -p          convert the base score to error probability\n"
		"   -n  <int>   read number for base content check [100000]\n"
		"   -e  <int>   read end length for base content check [150]\n"
		"   -5  <int>   trim bases from the 5' end of the read [auto]\n"
		"   -3  <int>   trim bases from the 3' end of the read [auto]\n"
		" Adapter filter options:\n"
		"   -a  <str>   adapter sequence file \n"
		"   -A          disable reads filter, only for adapter identify\n"
		"   -N  <int>   read number for adapter identify [100000]\n"
		"   -E  <int>   read end length for adapter trim [150]\n"
		"   -m  <int>   min match length for end adapter [4]\n"
		"   -M  <int>   min match length for middle adapter [35]\n"
		"   -T  <int>   extra trim length for middle adpter on both side [50]\n"
		"   -s <float>  min similarity for end adapter [0.75]\n"
		"   -S <float>  min similarity for middle adapter [0.9]\n"
		"   -D          discard reads with middle adapter instead of split\n"
		" Downsampling options:\n"
		"   -g  <str>   genome size (k/m/g)\n"
		"   -d  <int>   downsample to the desired coverage (requires -g) \n"
		"   -r  <int>   downsample to the desired number of reads \n"
		"   -R <float>  downsample to the desired fraction of reads \n"
		"   -F          disable reads filter, only for downsampling\n"
		" Other options:\n"
		"   -c  <int>   compression level (0-9) for compressed output [6]\n"
		"   -f          force FASTA output (discard quality) \n"
		"   -t  <int>   number of threads [16]\n"
		"   -h          show help [v1.09]\n"
		"\n";
	return 1;
}

int compLevel = 6;
int qType =0; //33 or 64 (old) for fastq

class Para_A24 {
	public:
		//
		string InFile;
		string OutFile;
		string TmpOutFile;
		//
		int MinLen;
		uint64_t MaxLen;
		int MinQ;
		uint64_t MaxQ;
		bool ErrorPro;
		int BCNum;
		int BCLen;
		int HeadCrop;
		int TailCrop;
		//
		string AdapterFile;
		int ReadNumber;
		int EndLen;
		int Kmer;
		int AdapterLen;
		int AdapterDep;
		int EndMatchLen;
		int MidMatchLen;
		int ExtraLen;
		float EndSim;
		float MidSim;
		bool discard;
		//
		uint64_t GenomeSize;
		int DesiredDepth;
		int DesiredNum;
		float DesiredFrac;
		bool Downsample;
		bool Filter;
		//
		bool FastaOut;
		int n_thread;
		int Infq;
		int Outfq;
		bool OUTGZ;
		bool ONLYAD;
		int ReadLength;

		Para_A24() {
			InFile="";
			OutFile="";
			
			MinLen=1000;
			MaxLen=UINT64_MAX;
			MinQ=0;
			MaxQ=UINT64_MAX;
			ErrorPro=false;
			BCNum=100000;
			BCLen=150;
			HeadCrop=-1;
			TailCrop=-1;
			
			AdapterFile="";
			ReadNumber=100000;
			EndLen=150;
			Kmer=19;
			AdapterLen=30;
			AdapterDep=100;
			EndMatchLen=4;
			MidMatchLen=35;
			ExtraLen=50;
			EndSim=0.75;
			MidSim=0.9;
			discard=false;

			GenomeSize=0;
			DesiredDepth=0;
			DesiredNum=0;
			DesiredFrac=0;
			Downsample=false;
			Filter=true;
			
			FastaOut=false;
			n_thread=16;
			Infq=3;
			Outfq=3;
			OUTGZ=false;
			ONLYAD=false;
			ReadLength=0;
		}
};

inline void error_exit(const string& msg) {
    cerr << "ERROR: " << msg << endl;
    exit(-1);
}

inline void  LogLackArg(string &flag) {
	cerr << "Error: Lack Argument for [ -"<<flag<<" ]"<<endl;
}

string & replace_all(string & str, const string & pattern, const string & replacement) {
	while(true) {
		string::size_type  pos(0);
		if((pos=str.find(pattern))!=string::npos){
			str.replace(pos,pattern.length(),replacement);
		} else{
			break;
		}
	}
	return str;
}

uint64_t GetGenomeSize(string &genomeSize) {
    std::string numberPart = genomeSize.substr(0, genomeSize.size() - 1);
    char suffix = genomeSize.back();
    double number = std::stod(numberPart);

    uint64_t result = 0;
    if (suffix == 'k' || suffix == 'K') {
        result = static_cast<uint64_t>(number * 1000.0);
    } else if (suffix == 'm' || suffix == 'M') {
        result = static_cast<uint64_t>(number * 1000000.0);
    } else if (suffix == 'g' || suffix == 'G') {
        result = static_cast<uint64_t>(number * 1000000000.0);
    } else {
        cerr << "Error: Genome size should end with k/m/g or K/M/G" << endl;
        return 0; 
    }

    return result;
}

int TGSFilter_cmd(int argc, char **argv, Para_A24 * P2In) {
	if (argc <= 2) {TGSFilter_usage(); return 1;}

	for(int i = 1; i < argc; i++){
		if(argv[i][0] != '-') {
			cerr << "Error: command option error! please check." << endl;
			return 1;
		}

		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		//input/output options
		if (flag == "i" ) {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->InFile=argv[i];
		}
		else if (flag == "o" ) {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->OutFile=argv[i];
		}

		//Basic filter options
		else if (flag == "l") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinLen=atoi(argv[i]);
			if (P2In->MinLen<100){
				P2In->MinLen=100;
			}
		}
		else if (flag == "L") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MaxLen=atoi(argv[i]);
		}
		else if (flag == "q") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinQ=atoi(argv[i]);
		}
		else if (flag == "Q") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MaxQ=atoi(argv[i]);
		}
		else if (flag == "p") {
			P2In->ErrorPro=true;
		}
		else if (flag == "n"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->BCNum=atoi(argv[i]);
		}
		else if (flag == "e"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->BCLen=atoi(argv[i]);
		}
		else if (flag == "5") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->HeadCrop=atoi(argv[i]);
		}
		else if (flag == "3"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->TailCrop=atoi(argv[i]);
		}
		
		//Adapter filter options
		else if (flag == "a"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->AdapterFile=argv[i];
		}
		else if (flag == "A"){
			P2In->ONLYAD=true;
		}
		else if (flag == "N"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->ReadNumber=atoi(argv[i]);
		}
		else if (flag == "E"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndLen=atoi(argv[i]);
		}
		else if (flag == "m"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndMatchLen=atoi(argv[i]);
		}
		else if (flag == "M"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MidMatchLen=atoi(argv[i]);
		}
		else if (flag == "T"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->ExtraLen=atoi(argv[i]);
		}
		else if (flag == "s"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndSim=atof(argv[i]);
			if (P2In->EndSim<0.7){
				P2In->EndSim=0.7;
				cerr << "Warning: re set -s to : "<<P2In->EndSim<<endl;
			}
		}
		else if (flag == "S"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MidSim=atof(argv[i]);
			if (P2In->MidSim<0.8){
				P2In->MidSim=0.8;
				cerr << "Warning: re set -S to : "<<P2In->MidSim<<endl;
			}
		}
		else if (flag == "D"){
			P2In->discard=true;
		}

		//Downsampling options
		else if (flag == "g"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			std::string genomeSize = argv[i];
			P2In->GenomeSize=GetGenomeSize(genomeSize);
			if(P2In->GenomeSize == 0){
				return 1;
			}
		}
		else if (flag == "d"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->DesiredDepth=atoi(argv[i]);
		}
		else if (flag == "r"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->DesiredNum=atoi(argv[i]);
		}
		else if (flag == "R"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->DesiredFrac=atof(argv[i]);
		}
		else if (flag == "F"){
			P2In->Filter=false;
		}

		//Other options
		else if (flag  ==  "c") {
			if(i + 1 == argc) {LogLackArg(flag) ; return 1;}
			i++;
			compLevel=atoi(argv[i]);
		}
		else if (flag  ==  "f") {
			P2In->FastaOut=true;
		}
		else if (flag  ==  "t") {
			if(i + 1 == argc) {LogLackArg(flag) ; return 1;}
			i++;
			P2In->n_thread=atoi(argv[i]);
		}
		else if (flag  == "help" || flag  == "h") {
			TGSFilter_usage(); return 1;
		}
		else {
			cerr << "Error: UnKnow argument -"<<flag<<endl;
			return 1;
		}
	}
	//check threads
	unsigned int maxThreads = std::thread::hardware_concurrency();
	if (maxThreads > 0 && P2In->n_thread > maxThreads){
		P2In->n_thread = maxThreads-2;
	}
	

	
	//check downsampling
	if (P2In->DesiredNum > 0 || P2In->DesiredFrac > 0){
		P2In->Downsample=true;
	}else{
		if (P2In->GenomeSize > 0 || P2In->DesiredDepth > 0){
			if (P2In->GenomeSize > 0 && P2In->DesiredDepth > 0){
				P2In->Downsample=true;
			}else if (P2In->GenomeSize > 0 && P2In->DesiredDepth == 0){
				cerr<< "The desired depth was required, along with the genome size!"<<endl;
				return 1;
			}else if (P2In->GenomeSize == 0 && P2In->DesiredDepth > 0){
				cerr<< "The genome size was required, along with the desired depth!"<<endl;
				return 1;
			}
		}else{
			P2In->Downsample=false;
		}
	}

	// check input and output
	if ((P2In->InFile).empty()) {
		cerr<< "Error: -i lack argument for the must"<<endl;
		return 1;
	}

	if (access((P2In->InFile).c_str(), 0) != 0) {
		cerr<<"Error: Can't find this file for -i "<<(P2In->InFile)<<endl;
		return 1;
	}
	
	return  0;
	
}

//read fasta or fastq file
#define FQ_BUF_SIZE (1<<21)
#define IGZIP_IN_BUF_SIZE (1<<20)
#define GZIP_HEADER_BYTES_REQ (1<<16)

struct ks {
    string name;
    string seq;
    string strand;
    string qual;
};

inline bool ends_with(string const & value,  string const & ending) {
	if (ending.size() > value.size()) return false;
	return  equal(ending.rbegin(), ending.rend(), value.rbegin());
}

class FastxReader {
public:
    FastxReader(string filename) {
        mFilename = filename;
        mZipped = false;
        isFastq = 2;
        mFile = NULL;
        mFastqBuf = new char[FQ_BUF_SIZE];
        mBufDataLen = 0;
        mBufUsedLen = 0;
        mGzipInputBufferSize = IGZIP_IN_BUF_SIZE;
        mGzipInputBuffer = new unsigned char[mGzipInputBufferSize];
        mGzipOutputBufferSize = FQ_BUF_SIZE;
        mGzipOutputBuffer = (unsigned char*)mFastqBuf;
        mGzipInputUsedBytes = 0;
        init();
    }

    ~FastxReader() {
        close();
        delete[] mFastqBuf;
        delete[] mGzipInputBuffer;
    }

    bool isZipped(){
        return mZipped;
    }

    ks* read(){
        if (isFastq == 1){
            return readFastq();
        }else{
            return readFasta();
        }
    }

private:
    void init(){
        mFile = fopen(mFilename.c_str(), "rb");
        if(mFile == NULL) {
            error_exit("Failed to open file: " + mFilename);
        }
        
        if (ends_with(mFilename, ".gz")){
            isal_gzip_header_init(&mGzipHeader);
            isal_inflate_init(&mGzipState);
            mGzipState.crc_flag = ISAL_GZIP_NO_HDR_VER;
            mGzipState.next_in = mGzipInputBuffer;
            mGzipState.avail_in = fread(mGzipState.next_in, 1, mGzipInputBufferSize, mFile);
            mGzipInputUsedBytes += mGzipState.avail_in;
            int ret = isal_read_gzip_header(&mGzipState, &mGzipHeader);
            if (ret != ISAL_DECOMP_OK) {
                error_exit("igzip: Error invalid gzip header found: "  + mFilename);
            }
            mZipped = true;
            if (ends_with(mFilename, ".fastq.gz") || ends_with(mFilename, ".fq.gz")){
                isFastq = 1;
            }else if (ends_with(mFilename, ".fasta.gz") || ends_with(mFilename, ".fa.gz")){
                isFastq = 0;
            }else{
                isFastq = 2;
            }
        }else{
            mZipped = false;
            if (ends_with(mFilename, ".fastq") || ends_with(mFilename, ".fq")){
                isFastq = 1;
            }else if (ends_with(mFilename, ".fasta") || ends_with(mFilename, ".fa")){
                isFastq = 0;
            }else{
                isFastq = 2;
            }
        }
        
        readToBuf();
    }

    void readToBuf(){
        mBufDataLen = 0;
        if(mZipped) {
            readToBufIgzip();
        } else {
            if(!eof())
                mBufDataLen = fread(mFastqBuf, 1, FQ_BUF_SIZE, mFile);
        }
        mBufUsedLen = 0;
    }

    void readToBufIgzip(){
        mBufDataLen = 0;
        while(mBufDataLen == 0) {
            if(eof() && mGzipState.avail_in==0)
                return;
            if (mGzipState.avail_in == 0) {
                mGzipState.next_in = mGzipInputBuffer;
                mGzipState.avail_in = fread(mGzipState.next_in, 1, mGzipInputBufferSize, mFile);
                mGzipInputUsedBytes += mGzipState.avail_in;
            }
            mGzipState.next_out = mGzipOutputBuffer;
            mGzipState.avail_out = mGzipOutputBufferSize;
            int ret = isal_inflate(&mGzipState);
            if (ret != ISAL_DECOMP_OK) {
                error_exit("igzip: encountered while decompressing file: " + mFilename);
            }
            mBufDataLen = mGzipState.next_out - mGzipOutputBuffer;
            if(eof() || mGzipState.avail_in>0)
                break;
        }
        // this block is finished
        if(mGzipState.block_state == ISAL_BLOCK_FINISH) {
            // a new block begins
            if(!eof() || mGzipState.avail_in > 0) {
                if (mGzipState.avail_in == 0) {
                    isal_inflate_reset(&mGzipState);
                    mGzipState.next_in = mGzipInputBuffer;
                    mGzipState.avail_in = fread(mGzipState.next_in, 1, mGzipInputBufferSize, mFile);
                    mGzipInputUsedBytes += mGzipState.avail_in;
                } else if (mGzipState.avail_in >= GZIP_HEADER_BYTES_REQ){
                    unsigned char* old_next_in = mGzipState.next_in;
                    size_t old_avail_in = mGzipState.avail_in;
                    isal_inflate_reset(&mGzipState);
                    mGzipState.avail_in = old_avail_in;
                    mGzipState.next_in = old_next_in;
                } else {
                    size_t old_avail_in = mGzipState.avail_in;
                    memmove(mGzipInputBuffer, mGzipState.next_in, mGzipState.avail_in);
                    size_t added = 0;
                    if(!eof()) {
                        added = fread(mGzipInputBuffer + mGzipState.avail_in, 1, mGzipInputBufferSize - mGzipState.avail_in, mFile);
                        mGzipInputUsedBytes += added;
                    }
                    isal_inflate_reset(&mGzipState);
                    mGzipState.next_in = mGzipInputBuffer;
                    mGzipState.avail_in = old_avail_in + added;
                }
                int ret = isal_read_gzip_header(&mGzipState, &mGzipHeader);
                if (ret != ISAL_DECOMP_OK) {
                    error_exit("igzip: invalid gzip header found");
                }
            }
        }

        if(eof() && mGzipState.avail_in == 0) {
            // all data was processed - fail if not at logical end of zip file (truncated?)
            if (mGzipState.block_state != ISAL_BLOCK_FINISH || !mGzipState.bfinal) {
                error_exit("igzip: unexpected eof");
            }
        }
    }

    bool bufferFinished(){
        if(mZipped) {
            return eof() && mGzipState.avail_in == 0;
        } else {
            return eof();
        }
    }

    void getLine(string* line){
        int copied = 0;
        int start = mBufUsedLen;
        int end = start;

        while(end < mBufDataLen) {
            if(mFastqBuf[end] != '\r' && mFastqBuf[end] != '\n'){
                end++;
            } else {
                break;
            }
        }

        // this line well contained in this buf, or this is the last buf
        if(end < mBufDataLen || bufferFinished()) {
            int len = end - start;
            line->assign(mFastqBuf+start, len);
            end++;
            // handle \r\n
            if(end < mBufDataLen-1 && mFastqBuf[end-1]=='\r' && mFastqBuf[end] == '\n'){
                end++;
            }
            mBufUsedLen = end;
            return ;
        }

        // this line is not contained in this buf, we need to read new buf
        line->assign(mFastqBuf+start, mBufDataLen - start);

        while(true) {
            readToBuf();
            start = 0;
            end = 0;
            while(end < mBufDataLen) {
                if(mFastqBuf[end] != '\r' && mFastqBuf[end] != '\n')
                    end++;
                else
                    break;
            }
            // this line well contained in this buf
            if(end < mBufDataLen || bufferFinished()) {
                int len = end - start;
                line->append(mFastqBuf+start, len);

                // skip \n or \r
                end++;
                // handle \r\n
                if(end < mBufDataLen-1 && mFastqBuf[end] == '\n')
                    end++;

                mBufUsedLen = end;
                return;
            }
            // even this new buf is not enough, although impossible
            line->append(mFastqBuf+start, mBufDataLen);
        }

        return;
    }

    ks* readFastq(){
        if (mBufUsedLen >= mBufDataLen && bufferFinished()) {
            return NULL;
        }

        ks* readObj = new ks();
        
        getLine(&(readObj->name));

        // name should start with @
        while ((readObj->name.empty() && !(mBufUsedLen >= mBufDataLen && bufferFinished())) 
                                || (!readObj->name.empty() && (readObj->name[0] != '@'))) { 
            getLine(&(readObj->name));
        }

		if (readObj->name.empty()) {
            delete readObj;
            return NULL;
        } else {
            readObj->name = readObj->name.substr(1);
        }

        getLine(&(readObj->seq));
        getLine(&(readObj->strand));
        getLine(&(readObj->qual));

        if (readObj->strand.empty() || (readObj->strand[0] != '+')) {
            std::cerr << "Error: third line should begin with '+', got: " << readObj->strand << endl;
            delete readObj;
            return NULL;
            
        }

        if (readObj->qual.empty() && readObj->seq.empty()) {
            cerr << "warning: sequence and quality are empty:" << readObj->name << endl;
            delete readObj;
            return NULL;
        }

        if (readObj->qual.length() != readObj->seq.length()) {
            cerr << "warning: sequence and quality have different length:" << readObj->name << endl;
            delete readObj;
            return NULL;
        }

        return readObj;
    }

    ks* readFasta(){
        if (mBufUsedLen >= mBufDataLen && bufferFinished()) {
            return NULL;
        }

        ks* readObj = new ks();
        
        getLine(&(readObj->name));

        // name should start with >
        while ((readObj->name.empty() && !(mBufUsedLen >= mBufDataLen && bufferFinished())) 
                                || (!readObj->name.empty() && (readObj->name[0] != '>'))) { 
            getLine(&(readObj->name));
        }

        if (readObj->name.empty()) {
            delete readObj;
            return NULL;
        } else {
            readObj->name = readObj->name.substr(1);
        }

        getLine(&(readObj->seq));

        if (readObj->seq.empty()) {
            cerr << "warning: sequence is empty:" << readObj->name << endl;
            delete readObj;
            return NULL;
        }

        return readObj;
    }

    bool eof(){
        return feof(mFile);
    }

    void close(){
        if (mFile){
		    fclose(mFile);
		    mFile = NULL;
	    }
    }

    string mFilename;
    struct isal_gzip_header mGzipHeader;
    struct inflate_state mGzipState;
    unsigned char *mGzipInputBuffer;
    unsigned char *mGzipOutputBuffer;
    size_t mGzipInputBufferSize;
    size_t mGzipOutputBufferSize;
    size_t mGzipInputUsedBytes;
    FILE* mFile;
    bool mZipped;
    int isFastq;
    char* mFastqBuf;
    int mBufDataLen;
    int mBufUsedLen;
};
//


class DeflateCompress {
public:
    DeflateCompress() {
        mCompressor = libdeflate_alloc_compressor(compLevel);
    }

    ~DeflateCompress() {
        libdeflate_free_compressor(mCompressor);
    }

    bool compressData(const void* input, uint8_t* &out, size_t & outSize) {
        size_t size = strlen(reinterpret_cast<const char*>(input));
        size_t bound = libdeflate_gzip_compress_bound(mCompressor, size);
        out = static_cast<uint8_t*>(malloc(bound));
        outSize = libdeflate_gzip_compress(mCompressor, input, size, out, bound);

        if (outSize == 0) {
            free(out);
            return false;
        } else {
            return true;
        }
    }

private:
    libdeflate_compressor* mCompressor;
};

string GetFileExtension(const string& FilePath) {
	size_t dotPos = FilePath.rfind('.');
	if (dotPos == string::npos) {
		return "";
	} 
	else {
		return FilePath.substr(dotPos + 1);
	}
}

string GetFilePreifx(const string& FilePath){
	string ext = GetFileExtension(FilePath);
	string prefix=FilePath;
	if (ext == "gz") {
		prefix=FilePath.substr(0, FilePath.rfind('.'));
		ext = GetFileExtension(prefix);
	}
	if (ext == "fq" || ext == "fastq" || ext == "fa" || ext == "fasta"){
		prefix=prefix.substr(0, prefix.rfind('.'));
	}
	return prefix;
}

int GetFileType(const string& FilePath) {
	int fqfile=3;
	string ext = GetFileExtension(FilePath);
	if (ext == "gz") {
		string FilePathA = FilePath.substr(0, FilePath.rfind('.'));
		ext = GetFileExtension(FilePathA);
	}

	if ((ext == "fa") || (ext == "fasta")) {
		fqfile = 0;
	} else if ((ext == "fq") || (ext == "fastq")) {
		fqfile = 1;
	}else if ((ext == "sam") || (ext == "SAM") || (ext == "bam") || (ext == "BAM")){
		fqfile = 2;
	}else{
		fqfile = 3;
	}
	return fqfile;
}

char complement[256];
string rev_comp_seq(const string& dna) {
	
	string reverse_complement;
	for (int i = dna.size() - 1; i >= 0; --i) {
		reverse_complement += complement[dna[i]];
	}
	return reverse_complement;
}

string adapterSearch(Para_A24 *P2In, std::vector<std::string>& reads, float &meanDep){
	/*
	PB-1 Pacific Biosciences Blunt Adapter
	PB-2  Pacific Biosciences C2 Primer
	ONT-1 Ligation
	ONT-2 Ligation
	ONT-3 Rapid
	ONT-4 1D^2
	ONT-5 1D^2
	ONT-6 Ligation(LA) / Native(NA) / Rapid(RA) / Rapid T(RAT) top strand
	ONT-7 Ligation Adapter bottom strand
	ONT-8 Native Adapter bottom strand
	ONT-9 cDNA RT Adapter (CRTA)
	*/
	std::map<std::string, std::pair<std::string, std::string>> adapterLib;

	adapterLib["PB-1"] = {"ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT", 
						 "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT"};
	adapterLib["PB-2"] = {"AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA", 
						 "TCCTCCTCCTCCGTTAATTTTTTTTTTTTTTTTTT"};
	adapterLib["ONT-1"] = {"AATGTACTTCGTTCAGTTACGTATTGCT", 
						  "AGCAATACGTAACTGAACGAAGTACATT"};
	adapterLib["ONT-2"] = {"GCAATACGTAACTGAACGAAGT",
						   "ACTTCGTTCAGTTACGTATTGC"};
	adapterLib["ONT-3"] = {"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA", 
						 	"TGAAGCGGCGCACGAAAAACGCGAAAGCGTTTCACGATAAATGCGAAAAC"};
	adapterLib["ONT-4"] = {"GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT",
						   "AGGAACCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCC"};
	adapterLib["ONT-5"] = {"GGAACCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCCAGCAAT", 
							"ATTGCTGGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCC"};
	adapterLib["ONT-6"] = {"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT",
						  "AGCAATACGTAACTGAACGAAGTACAGGAAAAAAAA"};
	adapterLib["ONT-7"] = {"GCAATACGTAACTGAACGAAGTACAGG",
						   "CCTGTACTTCGTTCAGTTACGTATTGC"};
	adapterLib["ONT-8"] = {"ACGTAACTGAACGAAGTACAGG", 
							"CCTGTACTTCGTTCAGTTACGT"};
	adapterLib["ONT-9"] = {"CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG",
						   "CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG"};
	//
	std::map<std::pair<std::string, std::string>, int> maps;

	float minSim=P2In->MidSim;
	for (const auto& ts : reads) {
        int tsLen=ts.length();
		string bestName;
		string bestStrand;
		int minLen=0;
        for (const auto& pair : adapterLib) {
            string name = pair.first;
			string qsFor=pair.second.first;
            string qsRev=pair.second.second;
            int qsForLen=qsFor.length();
            int qsRevLen=qsRev.length();
            int minK=static_cast<int>((1-minSim)*qsForLen)+1;
			//
            EdlibAlignResult result_for = edlibAlign(qsFor.c_str(), qsForLen, ts.c_str(), tsLen, 
                        edlibNewAlignConfig(minK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
            if (result_for.status == EDLIB_STATUS_OK){
                int dist=result_for.editDistance;
                int length=result_for.alignmentLength;
				int mlen=length-dist;
				int numAln=result_for.numLocations;
				if (numAln>0){
					if (mlen>minLen){
						bestName=name;
						bestStrand="+";
						minLen=mlen;
					}
				}
            }
            edlibFreeAlignResult(result_for);
            //
            EdlibAlignResult result_rev = edlibAlign(qsRev.c_str(), qsRevLen, ts.c_str(), tsLen, 
                        edlibNewAlignConfig(minK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
            if (result_rev.status == EDLIB_STATUS_OK){
                int dist=result_rev.editDistance;
                int length=result_rev.alignmentLength;
				int mlen=length-dist;
				int numAln=result_for.numLocations;
				if (numAln>0){
					if (mlen>minLen){
						bestName=name;
						bestStrand="-";
						minLen=mlen;
					}
				}
            }
            edlibFreeAlignResult(result_rev);
        }
		if (minLen>0){
			maps[std::make_pair(bestName, bestStrand)]+=minLen;
		}
    }

	//
	std::vector<std::pair<std::pair<std::string, std::string>, int>> maps_vec(maps.begin(), maps.end());
	std::sort(maps_vec.begin(), maps_vec.end(),
              [](const std::pair<std::pair<std::string, std::string>, int>& a,
                 const std::pair<std::pair<std::string, std::string>, int>& b) -> bool {
                  return a.second > b.second;
              });
	
	string name;
	string adapter;
	string strand;

	if (!maps_vec.empty()) {
        const auto& maxItem = maps_vec.front();
		name=maxItem.first.first;
		adapter=adapterLib[name].second;
		strand=maxItem.first.second;
		meanDep=static_cast<float>(maxItem.second)/(adapter.length());
		if (meanDep<2){
			adapter="";
			strand="";
			meanDep=0;
		}
    }

	if (strand=="-"){
		adapter=rev_comp_seq(adapter);
	}
	
	return adapter;
}

std::map<std::string, std::pair<std::string, std::string>> adapters;

int Get_filter_parameter(Para_A24 *P2In){
	string FilePath;
	string InPath=(P2In->InFile);
	string OutPath=(P2In->OutFile);
	if (OutPath.empty()){
		FilePath=InPath;
	}else{
		FilePath=OutPath;
	}

	int seqNum = 0;
	int ADLen = P2In->EndLen;
	if (ADLen > 100){
		ADLen = 100;
	}
	int BCLen = P2In->BCLen;
	int BCNum = P2In->BCNum;
	int ADNum = P2In->ReadNumber;
	int maxSeq = BCNum;
	if (maxSeq < ADNum){
		maxSeq = ADNum;
	}
	
	int qtNum=10000;
	int minQ=50000;
	int maxQ=0;

	std::vector<std::string> headSeq;
	std::vector<std::string> tailSeq;

	std::vector<int> headA(BCLen, 0);
	std::vector<int> headT(BCLen, 0);
	std::vector<int> tailA(BCLen, 0);
	std::vector<int> tailT(BCLen, 0);

	if ((P2In->Infq)==0 || (P2In->Infq)==1){
		FastxReader reader(InPath);
        ks* readObj = nullptr;
        while ((readObj = reader.read()) != nullptr) {
			if (seqNum>=maxSeq){
				break;
			}
			string name = readObj->name;
			string seq = readObj->seq;
			string qual = readObj->qual;
			delete readObj;

			int length=seq.length();


			if (length<5000){
				continue;
			}
			seqNum++;
			//
			if (seqNum < qtNum){
				int qual_len=qual.length();
				for (int i=0; i<qual_len; i++) {
					if(minQ>qual[i]) {
						minQ=qual[i];
					}
					if(maxQ<qual[i]) {
						maxQ=qual[i];
					}
				}
			}

			string head_seq=seq.substr(0, BCLen);
			string tail_seq=seq.substr(length-BCLen);
			if (seqNum<BCNum){
				for (int x=0; x < BCLen; x++){
					if (head_seq[x]=='A'){
						headA[x]++;
					}else if(head_seq[x]=='T'){
						headT[x]++;
					}

					if (tail_seq[x]=='A'){
						tailA[x]++;
					}else if(tail_seq[x]=='T'){
						tailT[x]++;
					}

				}
			}

			if (seqNum<ADNum){
				string head_ADSeq=seq.substr(0, ADLen);
				string tail_ADSeq=seq.substr(length-(ADLen));
				headSeq.push_back(head_ADSeq);
				tailSeq.push_back(tail_ADSeq);
			}
		}
	}else if ((P2In->Infq)=2){
		htsFile *bamFile = hts_open((P2In->InFile).c_str(), "r");
		if (bamFile == nullptr) {
			cerr << "Failed to open BAM file" <<P2In->InFile <<endl;
			return 1;
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();
		//
		while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			if (seqNum>=maxSeq){
				break;
			}
			int seqLen=bamRecord->core.l_qseq;
			if (seqLen<5000){
				continue;
			}
			seqNum++;

			uint8_t *seq = bam_get_seq(bamRecord);
			uint8_t *qual = bam_get_qual(bamRecord);

			string name = bam_get_qname(bamRecord);
			string sequence;
			string quality;

			for (int i = 0; i < seqLen; ++i) {
				char b=Base[bam_seqi(seq, i)];
				char q=qual[i] + 33;
				sequence+=b;
				quality+=q;
			}
			int qualLen=quality.length();
			//
			if (seqNum < qtNum){
				for (int i=0; i<qualLen; i++) {
					if(minQ>quality[i]) {
						minQ=quality[i];
					}
					if(maxQ<quality[i]) {
						maxQ=quality[i];
					}
				}
			}

			string head_seq=sequence.substr(0, BCLen);
			string tail_seq=sequence.substr(seqLen - BCLen);
			if (seqNum<BCNum){
				for (int x=0; x < BCLen; x++){
					if (head_seq[x]=='A'){
						headA[x]++;
					}else if(head_seq[x]=='T'){
						headT[x]++;
					}

					if (tail_seq[x]=='A'){
						tailA[x]++;
					}else if(tail_seq[x]=='T'){
						tailT[x]++;
					}

				}
			}

			if (seqNum<ADNum){
				string head_ADSeq=sequence.substr(0, ADLen);
				string tail_ADSeq=sequence.substr(seqLen-(ADLen));
				headSeq.push_back(head_ADSeq);
				tailSeq.push_back(tail_ADSeq);
			}
		}

		bam_destroy1(bamRecord);
    	bam_hdr_destroy(bamHeader);
    	hts_close(bamFile);
	}

	//base content
	std::vector<int> headBC;
	std::vector<int> tailBC;
	int maxBC=(P2In->ReadNumber)/100;
	for (int i = 0; i < BCLen; ++i) {
		if (abs(headA[i]-headT[i])>=maxBC){
			headBC.push_back(i+1);
		}

		if (abs(tailA[i]-tailT[i])>=maxBC){
			tailBC.insert(tailBC.begin(), BCLen-i);
		}
    }
	int headDrop=0;
	 if (!headBC.empty()){
		if(headBC.front()<=5){
			headDrop=headBC.front();
		}
		for (int x=0; x<headBC.size(); ++x){
			if ((headBC[x]-headDrop)<=5){
				headDrop=headBC[x];
			}
		}
	}

	int tailDrop=0;
	if (!tailBC.empty()){
		if(tailBC.front()<=5){
			tailDrop=tailBC.front();
		}
		for (int x=0; x<tailBC.size(); ++x){
			if ((tailBC[x]-tailDrop)<=5){
				tailDrop=tailBC[x];
			}
		}
	}

	if ((P2In->HeadCrop)<0){
		P2In->HeadCrop=headDrop;
	}
	if ((P2In->TailCrop)<0){
		P2In->TailCrop=tailDrop;
	}

	//mean phred quality
	if (maxQ>0){
		if(minQ >= 33 &&  minQ <= 78  &&  maxQ >= 33 && maxQ <= 127) {
			qType=33;
		}
		else if (minQ >= 64  &&  minQ <= 108  &&  maxQ >= 64 && maxQ <= 127){
			qType=64;
		}
		else if (minQ < 55) {
			qType=33;
		}
		else {
			qType=64;
		}
	}

	//search adapter
	if (!(P2In->AdapterFile).empty() && access((P2In->AdapterFile).c_str(), 0) == 0){
		string adapter_name=P2In->AdapterFile;
		FastxReader reader(adapter_name);
        ks* readObj = nullptr;
        while ((readObj = reader.read()) != nullptr) {
			string name = readObj->name;
			string seq_for = readObj->seq;
			string seq_rev=rev_comp_seq(seq_for);
			adapters[name] = {seq_for,seq_rev};
			cerr <<"INFO: adapter " << name <<" : "<<seq_for<< endl;
			delete readObj;
		}
	} else {
		string adapter_5p;
		string adapter_3p;
		float depth_5p=0;
		float depth_3p=0;
		int len_5p=0;
		int len_3p=0;
		std::vector<std::pair<std::string, int>> can5pAdapters;
		std::vector<std::pair<std::string, int>> can3pAdapters;

		cerr << "INFO: searching 5' adapter..."<<endl;
		adapter_5p=adapterSearch(P2In, headSeq, depth_5p);
		len_5p=adapter_5p.length();
		
		cerr << "INFO: searching 3' adapter..."<<endl;
		adapter_3p=adapterSearch(P2In, tailSeq, depth_3p);
		len_3p=adapter_3p.length();
		
		cerr <<"INFO: 5' adapter: " << adapter_5p << endl;
		cerr <<"INFO: 3' adapter: " << adapter_3p << endl;

		cerr <<"INFO: mean depth of 5' adapter: "<<depth_5p<<endl;
		cerr <<"INFO: mean depth of 3' adapter: "<<depth_3p<<endl;

		// store the adapter to adapters map
		len_5p=adapter_5p.length();
		len_3p=adapter_3p.length();

		if(len_5p != 0 || len_3p != 0){
			if (len_5p >0 && len_3p >0){
				string rev_3p=rev_comp_seq(adapter_3p);
				string rev_5p=rev_comp_seq(adapter_5p);
				if ((adapter_5p == adapter_3p) || (adapter_5p == rev_3p)){
					adapters["adapter_5-3p"]={adapter_5p, rev_5p};
				}else{
					adapters["adapter_5p"]={adapter_5p, rev_5p};
					adapters["adapter_3p"]={adapter_3p, rev_3p};
				}
			}else if(len_5p >0){
				string rev_5p=rev_comp_seq(adapter_5p);
				adapters["adapter_5p"]={adapter_5p, rev_5p};
			}else if(len_3p >0){
				string rev_3p=rev_comp_seq(adapter_3p);
				adapters["adapter_3p"]={adapter_3p, rev_3p};
			}
		}
	}

	//
	
	if (!(P2In->ONLYAD)){
		cerr << "INFO: trim 5' end length: "<<P2In->HeadCrop<<endl;
		cerr << "INFO: trim 3' end length: "<<P2In->TailCrop<<endl;
		cerr << "INFO: min output reads length: "<<P2In->MinLen<<endl;
		if ((P2In->Infq)==1 || (P2In->Infq)==2){
			cerr << "INFO: min Phred average quality score: "<<(P2In->MinQ)<<endl;
			//P2In->MinQ=(P2In->MinQ)+qType;
		}
	}
	
	return 0;
}

int GetEditDistance(Para_A24 *P2In, string &query, string &target, 
                    std::vector<std::vector<int>> &adapterRegions){
	
	int midNum=0;
	float endSim=P2In->EndSim;
	float midSim=P2In->MidSim;
	
	int qLen=query.length();
	int tLen=target.length();

	int end5Len=(P2In->EndLen) - (P2In->HeadCrop);
	int end3Len=(P2In->EndLen) - (P2In->TailCrop);
	if (end5Len<0){
		end5Len=0;
	}
	if (end3Len<0){
		end3Len=0;
	}
	int extraLen=(P2In->ExtraLen);

	int maxK=qLen-(P2In->MidMatchLen)+1;

	//check middle adapter
	int tsmLen=tLen- end5Len - end3Len;
	if (tsmLen>=qLen){
		string tsm = target.substr(end5Len, tsmLen);
		EdlibAlignResult result = edlibAlign(query.c_str(), qLen, tsm.c_str(), tsmLen, 
						edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result.status == EDLIB_STATUS_OK){
			int dist=result.editDistance;
			int numAln=result.numLocations;
			int length=result.alignmentLength;
			int mlen=length - dist;
			if (mlen >= (P2In->MidMatchLen)){
				for (int i=0; i<numAln; i++){
					int ts=result.startLocations[i] + end5Len;
					int te=result.endLocations[i] + end5Len + 1;
					float sim = static_cast<float>(mlen) / length;
					
					if (sim >= midSim){
						ts=ts - extraLen;
						te=te + extraLen;
						if (ts<0){ts=0;}
						if (te>tLen){te=tLen;}
						midNum++;
						adapterRegions.push_back({ts, te});
					}
				}
			}
		}
		edlibFreeAlignResult(result);
	}

	//check 5' and 3' end
	int checkLen= (P2In->EndLen) + int(qLen / endSim);
	if (checkLen > tLen){
		checkLen = tLen;
	}
	int ts5Len=checkLen - (P2In->HeadCrop);
	int ts3Len=checkLen - (P2In->TailCrop);

	//5' end
	if (ts5Len >= 5){
		maxK = ts5Len -5 + 1;
		string ts5 = target.substr(0, ts5Len);
		EdlibAlignResult result_ts5 = edlibAlign(query.c_str(), qLen, ts5.c_str(), ts5Len, 
					edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result_ts5.status == EDLIB_STATUS_OK){
			int dist=result_ts5.editDistance;
			int length=result_ts5.alignmentLength;
			int numAln=result_ts5.numLocations;
			int mlen=length - dist;
			if (mlen >= (P2In->EndMatchLen)){
				for (int i=0; i<numAln; i++){
					int ts=result_ts5.startLocations[i];
					int te=result_ts5.endLocations[i]+1;
					float sim1 = static_cast<float>(mlen) / te;
					float sim2 = static_cast<float>(mlen) / length;
					if (sim1 >= endSim || sim2 >= endSim){
						adapterRegions.push_back({0, te});
					}
				}
			}
		}
		edlibFreeAlignResult(result_ts5);
	}

	//3' end
	if (ts3Len >= 5){
		maxK = ts3Len -5 + 1;
		string ts3=target.substr(tLen-ts3Len);
		EdlibAlignResult result_ts3 = edlibAlign(query.c_str(), qLen, ts3.c_str(), ts3Len, 
					edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result_ts3.status == EDLIB_STATUS_OK){
			int dist=result_ts3.editDistance;
			int length=result_ts3.alignmentLength;
			int numAln=result_ts3.numLocations;
			int mlen=length - dist;
			if (mlen >= (P2In->EndMatchLen)){
				for (int i=0; i<numAln; i++){
					int ts=result_ts3.startLocations[i] + tLen - ts3Len;
					int te=result_ts3.endLocations[i]+1 + tLen - ts3Len;
					float sim1 = static_cast<float>(mlen) / (tLen-ts);
					float sim2 = static_cast<float>(mlen) / length;
					if (sim1 >= endSim || sim2 >= endSim){
						adapterRegions.push_back({ts, tLen});
					}
				}
				
			}
		}
		edlibFreeAlignResult(result_ts3);
	}
	return midNum;
}

void adapterMap(Para_A24 * P2In, string &rawSeq, int &rawLen,
				std::vector<std::vector<int>> &keepRegions,
                std::array<std::atomic<uint64_t>, 16> &DropInfo){

	std::vector<std::vector<int>> adapterRegions;
	int mid_for;
	int mid_rev;
	//
    for (const auto& pair : adapters) {
        string qsFor=pair.second.first;
        string qsRev=pair.second.second;
		mid_for = GetEditDistance(P2In, qsFor, rawSeq, adapterRegions);
		mid_rev = GetEditDistance(P2In, qsRev, rawSeq, adapterRegions);
    }

	int midNum=0;
	if (mid_for > 0 || mid_rev > 0){
		midNum=1;
	}
	DropInfo[14]+=midNum;

	if (midNum>0 && P2In->discard){
		DropInfo[15]+=rawLen;
	} else {
		if (adapterRegions.size() >= 2) {
			std::sort(adapterRegions.begin(), adapterRegions.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
				if (a[0] == b[0]) {
					return a[1] < b[1];
				}
				return a[0] < b[0];
			});
		}

		std::vector<std::vector<int>> mergedRegions;
		for (const auto& region : adapterRegions) {
			if (!mergedRegions.empty() && (mergedRegions.back()[1] >= region[0])) {
				mergedRegions.back()[1] = std::max(mergedRegions.back()[1], region[1]);
			} else {
				mergedRegions.push_back(region);
			}
		}

		int keep_len = 0;
		int currentStart = 0;
		if (mergedRegions.size()>=1){
			DropInfo[8]++;
			for (const auto& region : mergedRegions) {
				DropInfo[9]+=(region[1]-region[0]);
				if (region[0] > currentStart) {
					keep_len=region[0] - currentStart;
					if (keep_len >= (P2In->MinLen) && keep_len <= (P2In->MaxLen)){
						keepRegions.push_back({currentStart, keep_len});
					}else{
						DropInfo[12]++;
						DropInfo[13]+=keep_len;
					}
				}
				currentStart = region[1];
			}

			if (currentStart < rawLen) {
				keep_len=static_cast<int>(rawLen) - currentStart;
				if (keep_len >= (P2In->MinLen) && keep_len <= (P2In->MaxLen)){
					keepRegions.push_back({currentStart, keep_len});
				}else{
					DropInfo[12]++;
					DropInfo[13]+=keep_len;
				}
			}
		}else{
			if (rawLen >= (P2In->MinLen) && rawLen <= (P2In->MaxLen)){
				keepRegions.push_back({currentStart, rawLen});
			}else{
				DropInfo[12]++;
				DropInfo[13]+=rawLen;
			}
		}
	}
}

double CalcAvgQuality(const string &seq, const string &qual, 
                      std::vector<std::vector<uint64_t>> &baseQual,
                      std::vector<std::vector<uint64_t>> &baseCounts,
					  std::vector<int> &rawDiffQualBases) {
    int seqLen = seq.length();
    int qualLen = qual.length();
    if (seqLen == 0 || qualLen == 0 || seqLen != qualLen) {
        return 0.0;
    }

	if (baseQual.size() < seqLen || baseCounts.size() < seqLen){
		baseQual.resize(seqLen, std::vector<uint64_t>(5));
		baseCounts.resize(seqLen, std::vector<uint64_t>(5));
	}

	uint64_t sumQ = 0;
    int qValue = 0;
	char base;

    for (uint64_t i = 0; i < seqLen; i++) {
        qValue = qual[i] - qType;
        sumQ += qValue;
		base = seq[i];
		rawDiffQualBases[qValue]++;
		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
			baseQual[i][0]+=qValue;
		}else if (base== 'G' || base=='g'){
			baseCounts[i][1]++;
			baseQual[i][1]+=qValue;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][2]++;
			baseQual[i][2]+=qValue;
		} else if (base== 'T' || base=='t'){
			baseCounts[i][3]++;
			baseQual[i][3]+=qValue;
		}
		baseCounts[i][4]++;
		baseQual[i][4]+=qValue;
    }
    return static_cast<double>(sumQ) / seqLen;
}

void Get_tail_base_qual(const string &seq, const string &qual, 
                      std::vector<std::vector<uint64_t>> &baseQual,
                      std::vector<std::vector<uint64_t>> &baseCounts) {
    int seqLen = seq.length();
    int qualLen = qual.length();
    if (seqLen == 0 || qualLen == 0 || seqLen != qualLen) {
        return;
    }

	int length=150;
	if (length>seqLen) {
		length=seqLen;
	}

	if (baseQual.size() < length || baseCounts.size() < length){
		baseQual.resize(length, std::vector<uint64_t>(5));
		baseCounts.resize(length, std::vector<uint64_t>(5));
	}

	string tailSeq=seq.substr(seqLen-length);
	string tailQual=qual.substr(seqLen-length);

	int j = 0;
	int qValue = 0;
	char base;

    for (uint64_t i = 0; i < length; i++) {
		j=length-i-1;
		qValue = tailQual[j] - qType;
		base = tailSeq[j];

		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
			baseQual[i][0]+=qValue;
		}else if (base== 'G' || base=='g'){
			baseCounts[i][1]++;
			baseQual[i][1]+=qValue;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][2]++;
			baseQual[i][2]+=qValue;
		} else if (base== 'T' || base=='t'){
			baseCounts[i][3]++;
			baseQual[i][3]+=qValue;
		}
		baseCounts[i][4]++;
		baseQual[i][4]+=qValue;
    }
}

void Get_base_counts(const string &seq, 
                    std::vector<std::vector<uint64_t>> &baseCounts){
    int seqLen = seq.length();
    if (seqLen == 0){
        return;
    }

	if (baseCounts.size() < seqLen){
		baseCounts.resize(seqLen, std::vector<uint64_t>(5));
	}

    for (int i = 0; i < seqLen; i++){
		char base = seq[i];
		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
		}else if (base== 'G' || base=='g'){
			baseCounts[i][1]++;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][2]++;
		} else if (base== 'T' || base=='t'){
			baseCounts[i][3]++;
		}
		baseCounts[i][4]++;
    }
}

void Get_tail_base_counts(const string &seq, 
                    std::vector<std::vector<uint64_t>> &baseCounts){
    int seqLen = seq.length();
    if (seqLen == 0){
        return;
    }

	int length=150;
	if (length>seqLen) {
		length=seqLen;
	}
	if (baseCounts.size() < length){
		baseCounts.resize(length, std::vector<uint64_t>(5));
	}

	string tailSeq=seq.substr(seqLen-length);

	int j = 0;
	char base;

    for (int i = 0; i < length; i++){
		j=length-i-1;
		base = tailSeq[j];

		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
		}else if (base== 'G' || base=='g'){
			baseCounts[i][1]++;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][2]++;
		} else if (base== 'T' || base=='t'){
			baseCounts[i][3]++;
		}
		baseCounts[i][4]++;
    }
}

class TGSFilterTask {
public:
    TGSFilterTask(Para_A24 *P2In)
        : P2In(P2In), 
		  inputFastx_queue(), 
          inputBam_queue(), 
		  output_queue(), 
		  outputGz_queue(),
          headCrop(P2In->HeadCrop),
		  totalCrop((P2In->HeadCrop) + (P2In->TailCrop)),
          readNum(0), 
		  filterNum(0), 
		  read_done(false), 
		  filter_done(false),
		  outNum(0), 
		  writeNum(0), 
		  inQueueSize(0), 
		  outQueueSize(0), 
		  DropInfo(),
          seqLens(),
		  rawBases(0),
		  cleanBases(0), 
		  rawLens(), 
		  cleanLens(),
		  worker_count(P2In->n_thread),
		  rawDiffQualBases(P2In->n_thread, std::vector<int>(256,0)),
		  rawDiffQualReads(P2In->n_thread, std::vector<int>(256,0)),
		  cleanDiffQualBases(P2In->n_thread, std::vector<int>(256,0)),
		  cleanDiffQualReads(P2In->n_thread, std::vector<int>(256,0)),
          rawBaseQual(P2In->n_thread),
		  rawBaseCounts(P2In->n_thread),
		  raw3pBaseQual(P2In->n_thread),
		  raw3pBaseCounts(P2In->n_thread),
	  	  cleanBaseQual(P2In->n_thread),
	  	  cleanBaseCounts(P2In->n_thread),
		  clean3pBaseQual(P2In->n_thread),
		  clean3pBaseCounts(P2In->n_thread) {}

    std::array<std::atomic<uint64_t>, 16> DropInfo;
    std::unordered_map<std::string, int> seqLens;

	uint64_t rawBases, cleanBases;
	std::vector<int> rawLens, cleanLens;

	std::vector<std::vector<int>> rawDiffQualBases, rawDiffQualReads;
	std::vector<std::vector<int>> cleanDiffQualBases, cleanDiffQualReads;
	
	std::vector<std::vector<std::vector<uint64_t>>> rawBaseQual, rawBaseCounts, raw3pBaseQual, raw3pBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> cleanBaseQual, cleanBaseCounts, clean3pBaseQual, clean3pBaseCounts;

    void start() {
		std::vector<std::thread> read_threads; // read
		if ((P2In->Infq)==2){
			read_threads.emplace_back(&TGSFilterTask::read_bam, this);
		} else if((P2In->Infq)==1 || (P2In->Infq)==0){
			read_threads.emplace_back(&TGSFilterTask::read_fastx, this);
		}

		std::vector<std::thread> filter_threads; // filter
		for (int i = 0; i < worker_count; ++i) {
			filter_threads.emplace_back(&TGSFilterTask::filter_sequence, this, i);
		}

		//
		std::vector<std::thread> output_threads; // output
		if (P2In->OUTGZ && !(P2In->Downsample)){
			output_threads.emplace_back(&TGSFilterTask::write_output_gz, this);
		} else {
			output_threads.emplace_back(&TGSFilterTask::write_output, this);
		}

		// Wait for all threads to finish
		for (auto& read_thread : read_threads) {
			read_thread.join();
		}
		for (auto& filter_thread : filter_threads) {
			filter_thread.join();
		}
		for (auto& output_thread : output_threads) {
			output_thread.join();
		}
	}


private:
    int read_fastx() {
		//
		string InPath = P2In->InFile;
		FastxReader reader(InPath);
        ks* readObj = nullptr;
        while ((readObj = reader.read()) != nullptr) {
			string name = readObj->name;
			string seq = readObj->seq;
			string qual = readObj->qual;
			int seqLen = seq.length();
			rawBases += seqLen;
			delete readObj;
			
			DropInfo[0]++;
			DropInfo[1] += seqLen;

			rawLens.push_back(seqLen);

			int checkLen=seqLen - totalCrop;
            if (checkLen >= (P2In->MinLen) && checkLen <= (P2In->MaxLen)){
				inputFastx_queue.enqueue(std::make_tuple(name, seq, qual));
				inQueueSize++;
				readNum++;
				if (inQueueSize >= 2 * worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
            } else {
				DropInfo[4]++;
				DropInfo[5] += seqLen ;
			}  
        }
        read_done = true;
        return 0;
    }

    int read_bam() {
        htsFile *bamFile = hts_open((P2In->InFile).c_str(), "r");
		if (bamFile == nullptr) {
			std::cerr << "Error: Failed to open file: " <<P2In->InFile <<endl;
			return 1;
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			string name = bam_get_qname(bamRecord);
			uint8_t *seq = bam_get_seq(bamRecord);
			uint8_t *qual = bam_get_qual(bamRecord);

			int seqLen = bamRecord->core.l_qseq;
			rawBases += seqLen;
			rawLens.push_back(seqLen);

			int checkLen=seqLen - totalCrop;

			if (checkLen >= (P2In->MinLen) && checkLen <= (P2In->MaxLen)){
				inputBam_queue.enqueue(std::make_tuple(name, seq, qual, seqLen));
				inQueueSize++;
				readNum++;
				if (inQueueSize >= 2 * worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}else{
				DropInfo[4]++;
				DropInfo[5] += seqLen ;
			}
		}
        read_done = true;
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
        return 0;
    }

    //
    void filter_sequence(int tid) {
        while (!(read_done && inQueueSize == 0 && readNum == filterNum)) {
			string rawName;
            string rawSeq;
            string rawQual;
			int rawSeqLen = 0;
			int rawQualLen=0;
			if ((P2In->Infq)==2){
				std::tuple<std::string, uint8_t*, uint8_t*, int> info;
				if (inputBam_queue.try_dequeue(info)) {
					inQueueSize--;
					rawName = std::get<0>(info);
					uint8_t* seqChar = std::get<1>(info);
					uint8_t* qualChar = std::get<2>(info);
					rawSeqLen = std::get<3>(info);
					rawQualLen = rawSeqLen;
					for (int i = 0; i < rawSeqLen; ++i) {
						char b = Base[bam_seqi(&seqChar[i], 0)];
						char q = qualChar[i] + 33;
						rawSeq += b;
						rawQual += q;
					}
				} else {
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			} else if ((P2In->Infq)==1 || (P2In->Infq)==0){
				std::tuple<std::string, std::string, std::string> info;
				if (inputFastx_queue.try_dequeue(info)) {
					inQueueSize--;
					rawName = std::get<0>(info);
					rawSeq = std::get<1>(info);
					rawQual = std::get<2>(info);
					rawSeqLen = rawSeq.length();
					rawQualLen = rawQual.length();
				} else {
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}
			
			if (rawSeqLen > 0){
				if (rawQualLen > 0){
					double avgQuality;
					avgQuality = CalcAvgQuality(rawSeq, rawQual, rawBaseQual[tid], rawBaseCounts[tid], rawDiffQualBases[tid]);
					rawDiffQualReads[tid][int(avgQuality)]++;
					Get_tail_base_qual(rawSeq, rawQual, raw3pBaseQual[tid], raw3pBaseCounts[tid]);
					if (avgQuality < (P2In->MinQ)){
						DropInfo[10]++;
						DropInfo[11]+= rawSeqLen;
						continue;
					}
					rawSeq = rawSeq.substr(headCrop, rawSeqLen - totalCrop);
                    rawQual = rawQual.substr(headCrop, rawQualLen - totalCrop);
				    rawSeqLen = rawSeq.length();
				    rawQualLen = rawQual.length();
				} else {
                    Get_base_counts(rawSeq, rawBaseCounts[tid]);
					Get_tail_base_counts(rawSeq, raw3pBaseCounts[tid]);
				    rawSeq = rawSeq.substr(headCrop, rawSeqLen - totalCrop);
				    rawSeqLen = rawSeq.length();
                }

				std::vector<std::vector<int>> keepRegions;
				adapterMap(P2In, rawSeq, rawSeqLen, keepRegions, DropInfo);

				string name;
				string seq;
                string qual;
				string out;
				
				int passNum=0;

				if (keepRegions.size()>0){
					for (const auto& region : keepRegions) {
						passNum++;
						int start = region[0];
						int seqLen = region[1];
						seq = rawSeq.substr(start,seqLen);

                        if (rawQualLen > 0){
                            qual=rawQual.substr(start,seqLen);
							double avgQuality;
                            avgQuality = CalcAvgQuality(seq, qual, cleanBaseQual[tid], cleanBaseCounts[tid], cleanDiffQualBases[tid]);
							cleanDiffQualReads[tid][int(avgQuality)]++;
							Get_tail_base_qual(seq, qual, clean3pBaseQual[tid], clean3pBaseCounts[tid]);
                        }else{
                            Get_base_counts(seq, cleanBaseCounts[tid]);
							Get_tail_base_counts(seq, clean3pBaseCounts[tid]);
                        }

                        //
						DropInfo[2]++;
						DropInfo[3]+=seqLen;

						outNum++;

						if (passNum>=2){
							name = rawName + ":" + std::to_string(passNum);
						}else{
							name=rawName;
						}

						if ((P2In->Outfq)==1){
							out = "@" + name +"\n" + seq + "\n+\n" + qual + "\n";
							if (P2In->OUTGZ && !(P2In->Downsample)){
								uint8_t *ComData;
								size_t ComSize;
								DeflateCompress GZData;
								if (GZData.compressData(out.c_str(), ComData, ComSize)) {
									outputGz_queue.enqueue(std::make_tuple(ComData, ComSize, seqLen));
									outQueueSize++;
								} else {
									free(ComData);
								}
							} else{
								output_queue.enqueue(std::make_tuple(out, name, seqLen));
								outQueueSize++;

							}
						} else if((P2In->Outfq)==0){
							out = ">" + name +"\n" + seq + "\n";
							if (P2In->OUTGZ && !(P2In->Downsample)){
								uint8_t *ComData;
								size_t ComSize;
								DeflateCompress GZData;
								if (GZData.compressData(out.c_str(), ComData, ComSize)) {
									outputGz_queue.enqueue(std::make_tuple(ComData, ComSize, seqLen));
									outQueueSize++;
								} else {
									free(ComData);
								}
							} else {
								output_queue.enqueue(std::make_tuple(out, name, seqLen));
								outQueueSize++;
							}
						}

						if (outQueueSize >= 2 * worker_count){
							this_thread::sleep_for(chrono::milliseconds(1));
						}
					}
				}
				filterNum++;
			} 
        }
		filter_done = true;
    }
    //

    void write_output_gz() {
        std::ofstream OUTHGZ((P2In->OutFile).c_str(), std::ios::out | std::ios::binary);
        if (!OUTHGZ.is_open()) {
			std::cerr << "Error: Failed to open file: " <<P2In->OutFile <<endl;
			return ;
        }

        while (!(read_done && filter_done && inQueueSize == 0 && outNum == writeNum && outQueueSize == 0 && readNum == filterNum )) {
			std::tuple<uint8_t*, size_t, int> info;
            if (outputGz_queue.try_dequeue(info)) {
				uint8_t* ComData=std::get<0>(info);
				size_t ComSize=std::get<1>(info);
                int seqLen=std::get<2>(info);
				cleanBases += seqLen;
				cleanLens.push_back(seqLen);
				if (ComData != nullptr && ComSize > 0) {
					OUTHGZ.write(reinterpret_cast<const char*>(ComData), ComSize);
					free(ComData);
					writeNum++;
					outQueueSize--;
				}
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
        OUTHGZ.close();
    }
    //
    void write_output() {
		if ((P2In->TmpOutFile).empty() && (P2In->OutFile).empty()){
			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
				std::tuple<std::string, std::string, int> info;
				if (output_queue.try_dequeue(info)) {
					string out=std::get<0>(info);
					string name=std::get<1>(info);
					int seqLen=std::get<2>(info);
					cleanBases += seqLen;
					cleanLens.push_back(seqLen);
					seqLens[name] = seqLen;
					std::cout << out;
					writeNum++;
					outQueueSize--;
				} else {
					std::this_thread::sleep_for(std::chrono::milliseconds(1));
				}
			}
		}else{
			string outName;
			if ((P2In->TmpOutFile).empty()){
				outName = (P2In->OutFile);
			}else{
				outName = (P2In->TmpOutFile);
			}

			std::ofstream OUTH(outName.c_str());
			if (!OUTH.is_open()) {
				std::cerr << "Error: Failed to open file: " << outName <<endl;
				return ;
			}

			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
				std::tuple<std::string, std::string, int> info;
				if (output_queue.try_dequeue(info)) {
					string out=std::get<0>(info);
					string name=std::get<1>(info);
					int seqLen=std::get<2>(info);
					cleanBases += seqLen;
					cleanLens.push_back(seqLen);
					seqLens[name] = seqLen;
					OUTH << out;
					writeNum++;
					outQueueSize--; 
				} else {
					std::this_thread::sleep_for(std::chrono::milliseconds(1));
				}
			}
			OUTH.close();
		}
	}

	void write_stdout() {
        
    }

    Para_A24 *P2In;
    moodycamel::ConcurrentQueue<std::tuple<std::string, std::string, std::string>> inputFastx_queue;
    moodycamel::ConcurrentQueue<std::tuple<std::string, uint8_t*, uint8_t*, int>> inputBam_queue;
	moodycamel::ConcurrentQueue<std::tuple<std::string, std::string, int>> output_queue;
    moodycamel::ConcurrentQueue<std::tuple<uint8_t*, size_t, int>> outputGz_queue;

    int headCrop;
    int totalCrop;

    std::atomic<bool> read_done;
	std::atomic<bool> filter_done;
	std::atomic<int> readNum;
	std::atomic<int> filterNum;
	std::atomic<int> outNum;
	std::atomic<int> writeNum;
	std::atomic<int> inQueueSize;
	std::atomic<int> outQueueSize;

	int worker_count;
};

class DownSampleTask {
public:
    DownSampleTask(Para_A24 *P2In, string &input, std::unordered_map<std::string, int> &seqLens)
        : P2In(P2In), 
		  InPath(input), 
		  seqLens(seqLens), 
		  worker_count(P2In->n_thread), 
		  totalSize(0),
		  seqNames(), 
		  inputFastx_queue(), 
		  inputBam_queue(), 
		  output_queue(), 
		  outputGz_queue(),
          OutInfo(), 
		  outNum(0), 
		  writeNum(0), 
		  readNum(0), 
		  filterNum(0), 
		  read_done(false), 
		  filter_done(false),
		  inQueueSize(0), 
		  outQueueSize(0),
		  downLens(),
		  downBases(0),
		  downDiffQualBases(P2In->n_thread, std::vector<int>(256,0)),
		  downDiffQualReads(P2In->n_thread, std::vector<int>(256,0)),
		  downBaseQual(P2In->n_thread),
	  	  downBaseCounts(P2In->n_thread),
		  down3pBaseQual(P2In->n_thread),
		  down3pBaseCounts(P2In->n_thread) {}
	
	std::array<std::atomic<uint64_t>, 2> OutInfo;

	uint64_t downBases;
	std::vector<int> downLens;
	std::vector<std::vector<int>> downDiffQualBases, downDiffQualReads;
	std::vector<std::vector<std::vector<uint64_t>>> downBaseQual, downBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> down3pBaseQual, down3pBaseCounts;


    void start() {
		//get output names
		if (!(P2In->Filter)){
			if ((P2In->Infq)==2){
				get_bam_SeqLen();
			}else if((P2In->Infq)==1 || (P2In->Infq)==0){
				get_fastx_SeqLen();
			}
		}
		get_reads_name();
		//
		std::vector<std::thread> read_threads; // read
		if (!(P2In->Filter)){
			if ((P2In->Infq)==2){
				read_threads.emplace_back(&DownSampleTask::read_bam, this);
			} else if((P2In->Infq)==1 || (P2In->Infq)==0){
				read_threads.emplace_back(&DownSampleTask::read_fastx, this);
			}
		}else{
			read_threads.emplace_back(&DownSampleTask::read_fastx, this);
		}

		std::vector<std::thread> filter_threads; // filter
		for (int i = 0; i < worker_count; ++i) {
			filter_threads.emplace_back(&DownSampleTask::filter_sequence, this, i);
		}

		//
		std::vector<std::thread> output_threads; // output
		if (P2In->OUTGZ){
			output_threads.emplace_back(&DownSampleTask::write_output_gz, this);
		} else {
			output_threads.emplace_back(&DownSampleTask::write_output, this);
		}

		// Wait for all threads to finish
		
		for (auto& read_thread : read_threads) {
			read_thread.join();
		}
		for (auto& filter_thread : filter_threads) {
			filter_thread.join();
		}
		for (auto& output_thread : output_threads) {
			output_thread.join();
		}
	}


private:
	void get_fastx_SeqLen(){
		FastxReader reader(InPath);
        ks* readObj = nullptr;
        while ((readObj = reader.read()) != nullptr) {
			string name = readObj->name;
			string seq = readObj->seq;
			int seqLen = seq.length();
			totalSize += seqLen;
			delete readObj;
			seqLens[name] = seqLen;
		}
	}

	int get_bam_SeqLen(){
		htsFile *bamFile = hts_open(InPath.c_str(), "r");
		if (bamFile == nullptr) {
			error_exit("Failed to open file: " + InPath);
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			string name = bam_get_qname(bamRecord);
			int seqLen = bamRecord->core.l_qseq;
			seqLens[name] = seqLen;
			totalSize += seqLen;
		}
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
		return 0;
	}

	void get_reads_name(){
		std::vector<std::pair<std::string, int>> vec(seqLens.begin(), seqLens.end());
		std::sort(vec.begin(), vec.end(), [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
    		return a.second > b.second;
		});

		if ((P2In->GenomeSize)>0 && (P2In->DesiredDepth)>0){
			uint64_t addedSize=0;
			uint64_t desiredSize=(P2In->GenomeSize) * (P2In->DesiredDepth);
			for(const auto& pair : vec) {
				seqNames.insert(pair.first);
				addedSize += pair.second;
				OutInfo[0]++;
				OutInfo[1]+=pair.second;
				downBases += pair.second;
				downLens.push_back(pair.second);
				if (addedSize >= desiredSize){
					break;
				}
			}
		}else if ((P2In->DesiredFrac)>0){
			if (totalSize==0){
				for(const auto& pair : vec) {
					totalSize += pair.second;
				}
			}
			uint64_t addedSize=0;
			uint64_t desiredSize=(P2In->DesiredFrac) * totalSize;
			for(const auto& pair : vec) {
				seqNames.insert(pair.first);
				addedSize += pair.second;
				OutInfo[0]++;
				OutInfo[1]+=pair.second;
				downBases += pair.second;
				downLens.push_back(pair.second);
				if (addedSize >= desiredSize){
					break;
				}
			}
		}else if ((P2In->DesiredNum)>0){
			int addedNum=0;
			for(const auto& pair : vec) {
				seqNames.insert(pair.first);
				addedNum++;
				OutInfo[0]++;
				OutInfo[1]+=pair.second;
				downBases += pair.second;
				downLens.push_back(pair.second);
				if (addedNum >= (P2In->DesiredNum)){
					break;
				}
			}
		}
	}

    int read_fastx() {
		FastxReader reader(InPath);
        ks* readObj = nullptr;

        while ((readObj = reader.read()) != nullptr) {
			string name = readObj->name;
			string seq = readObj->seq;
			string qual = readObj->qual;
			delete readObj;

			if (seqNames.find(name) != seqNames.end()) {
				inputFastx_queue.enqueue(std::make_tuple(name, seq, qual));
				inQueueSize++;
				readNum++;
				if (inQueueSize >= 2 * worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}
        }
		read_done = true;

        return 0;
    }

    int read_bam() {
        htsFile *bamFile = hts_open((P2In->InFile).c_str(), "r");
		if (bamFile == nullptr) {
			std::cerr << "Error: Failed to open file: " <<P2In->InFile <<endl;
			return 1;
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			string name = bam_get_qname(bamRecord);
			uint8_t *seq = bam_get_seq(bamRecord);
			uint8_t *qual = bam_get_qual(bamRecord);
			int seqLen = bamRecord->core.l_qseq;

			if (seqNames.find(name) != seqNames.end()) {
				inputBam_queue.enqueue(std::make_tuple(name, seq, qual, seqLen));
				inQueueSize++;
				readNum++;
				if (inQueueSize >= worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}
		}
        read_done = true;
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
        return 0;
    }

    //
    void filter_sequence(int tid) {
        while (!(read_done && inQueueSize == 0 && readNum == filterNum)) {
			string rawName;
            string rawSeq;
            string rawQual;
            int rawQualLen=0;
			int rawSeqLen = 0;

			if ((P2In->Infq)==2){
				std::tuple<std::string, uint8_t*, uint8_t*, int> info;
				if (inputBam_queue.try_dequeue(info)) {
					inQueueSize--;
					outNum++;
					rawName = std::get<0>(info);
					uint8_t* seqChar = std::get<1>(info);
					uint8_t* qualChar = std::get<2>(info);
					rawSeqLen = std::get<3>(info);
					rawQualLen = rawSeqLen;
					for (int i = 0; i < rawSeqLen; ++i) {
						char b = Base[bam_seqi(&seqChar[i], 0)];
						char q = qualChar[i] + 33;
						rawSeq += b;
						rawQual += q;
					}
				} else {
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			} else if ((P2In->Infq)==1 || (P2In->Infq)==0){
				std::tuple<std::string, std::string, std::string> info;
				if (inputFastx_queue.try_dequeue(info)) {
					inQueueSize--;
					outNum++;
					rawName = std::get<0>(info);
					rawSeq = std::get<1>(info);
					rawQual = std::get<2>(info);
					rawSeqLen = rawSeq.length();
					rawQualLen = rawQual.length();
				} else {
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}

			if (rawSeqLen > 0){
				if (rawQualLen > 0){
					double avgQuality;
                    avgQuality=CalcAvgQuality(rawSeq, rawQual, downBaseQual[tid], downBaseCounts[tid], downDiffQualBases[tid]);
					downDiffQualReads[tid][int(avgQuality)]++;
					Get_tail_base_qual(rawSeq, rawQual, down3pBaseQual[tid], down3pBaseCounts[tid]);
                }else{
                    Get_base_counts(rawSeq, downBaseCounts[tid]);
					Get_tail_base_counts(rawSeq, down3pBaseCounts[tid]);
                }

				string out;
				if ((P2In->Outfq)==1){
					out = "@" + rawName +"\n" + rawSeq + "\n+\n" + rawQual + "\n";
					if (P2In->OUTGZ){
						uint8_t *ComData;
						size_t ComSize;
						DeflateCompress GZData;
						if (GZData.compressData(out.c_str(), ComData, ComSize)) {
							outputGz_queue.enqueue({ComData, ComSize});
							outQueueSize++;
						} else {
							free(ComData);
						}
					} else{
						output_queue.enqueue(out);
						outQueueSize++;

					}
				} else if((P2In->Outfq)==0){
					out = ">" + rawName +"\n" + rawSeq + "\n";
					if (P2In->OUTGZ){
						uint8_t *ComData;
						size_t ComSize;
						DeflateCompress GZData;
						if (GZData.compressData(out.c_str(), ComData, ComSize)) {
							outputGz_queue.enqueue({ComData, ComSize});
							outQueueSize++;
						} else {
							free(ComData);
						}
					} else {
						output_queue.enqueue(out);
						outQueueSize++;
					}
				}

				if (outQueueSize >= worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
				filterNum++;
			}
        }
		filter_done = true;
    }
    //

    void write_output_gz() {
        std::ofstream OUTHGZ((P2In->OutFile).c_str(), std::ios::out | std::ios::binary);
        if (!OUTHGZ.is_open()) {
            throw std::runtime_error("Failed to open output file");
        }

        while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
            std::pair<uint8_t*, size_t> out;
            if (outputGz_queue.try_dequeue(out)) {
                if (out.first != nullptr && out.second > 0) {
                    OUTHGZ.write(reinterpret_cast<const char*>(out.first), out.second);
                    free(out.first);
					writeNum++;
					outQueueSize--;
                }
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }

        OUTHGZ.close();
    }
    //
    void write_output() {
		if (!((P2In->OutFile).empty())){
			std::ofstream OUTH((P2In->OutFile).c_str());
        	if (!OUTH.is_open()) {
            	throw std::runtime_error("Failed to open output file");
        	}
			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
            	std::string out;
            	if (output_queue.try_dequeue(out)) {
                	OUTH << out;
					writeNum++;
					outQueueSize--;
            	} else {
                	std::this_thread::sleep_for(std::chrono::milliseconds(1));
            	}
       	 	}

        	OUTH.close();
		}else{
			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
            	std::string out;
            	if (output_queue.try_dequeue(out)) {
					std::cout << out;
					writeNum++;
					outQueueSize--;
            	} else {
                	std::this_thread::sleep_for(std::chrono::milliseconds(1));
            	}
        	}
		}
	}

    Para_A24 *P2In;
	string InPath;
	std::unordered_map<std::string, int> seqLens;
	uint64_t totalSize;
	std::unordered_set<std::string> seqNames;
    int worker_count;
    moodycamel::ConcurrentQueue<std::tuple<std::string, std::string, std::string>> inputFastx_queue;
    moodycamel::ConcurrentQueue<std::tuple<std::string, uint8_t*, uint8_t*, int>> inputBam_queue;
    moodycamel::ConcurrentQueue<std::string> output_queue;
    moodycamel::ConcurrentQueue<std::pair<uint8_t*, int>> outputGz_queue;

    int headCrop;
    int totalCrop;

	std::atomic<bool> read_done;
	std::atomic<bool> filter_done;
	std::atomic<int> readNum;
	std::atomic<int> filterNum;
	std::atomic<int> outNum;
	std::atomic<int> writeNum;
	std::atomic<int> inQueueSize;
	std::atomic<int> outQueueSize;
};


std::string generateRandomString(const int length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis('a', 'z');

    std::string randomString;
    for (int i = 0; i < length; ++i) {
        randomString += static_cast<char>(dis(gen));
    }

    return randomString;
}

struct BarPlotData {
    std::vector<int> x;
    std::vector<float> y;
};

struct LineYData {
    std::string bases;
    std::vector<float> data;
};

struct LinePlotData {
    std::vector<int> x;
    std::vector<LineYData> y;
};

void Get_diff_quals(const std::vector<std::vector<int>> &diffQuals,
                    std::vector<uint64_t> &outQuals) {
    uint64_t Q5 = 0;
    uint64_t Q7 = 0;
    uint64_t Q10 = 0;
    uint64_t Q13 = 0;
    uint64_t Q15 = 0;
    uint64_t Q17 = 0;
    uint64_t Q20 = 0;
    uint64_t Q23 = 0;
    uint64_t Q30 = 0;
	uint64_t Q40 = 0;

    for (int i = 0; i < diffQuals.size(); i++) {
        for (int j = 0; j < diffQuals[i].size(); j++) {
            uint64_t number = diffQuals[i][j];
            if (j >= 5) {
                Q5 += number;
            }
            if (j >= 7) {
                Q7 += number;
            }
            if (j >= 10) {
                Q10 += number;
            }
            if (j >= 13) {
                Q13 += number;
            }
            if (j >= 15) {
                Q15 += number;
            }
            if (j >= 17) {
                Q17 += number;
            }
            if (j >= 20) {
                Q20 += number;
            }
            if (j >= 23) {
                Q23 += number;
            }
            if (j >= 30) {
                Q30 += number;
            }
			if (j >= 40) {
                Q40 += number;
            }
        }
    }
    outQuals = {Q5, Q7, Q10, Q13, Q15, Q17, Q20, Q23, Q30, Q40};
}

void Get_plot_line_data(Para_A24* P2In,
                        const int &rawMax, 
                        float &rawGC,
                        float &rawMeanQual,
                        const uint64_t &rawBases,
                        LinePlotData &rawReadsQual,
                        LinePlotData &rawBasesContents,
						LinePlotData &raw5pReadsQual,
                        LinePlotData &raw5pBasesContents,
						LinePlotData &raw3pReadsQual,
                        LinePlotData &raw3pBasesContents,
                        const std::vector<std::vector<std::vector<uint64_t>>> &rawBaseQual,
                        const std::vector<std::vector<std::vector<uint64_t>>> &rawBaseCounts,
						const std::vector<std::vector<std::vector<uint64_t>>> &raw3pBaseQual,
                        const std::vector<std::vector<std::vector<uint64_t>>> &raw3pBaseCounts){

    int n_threads=P2In->n_thread;
    uint64_t rawBaseQualSum=0;
    uint64_t rawGCSum=0;
    std::vector<std::string> bases = {"A", "G", "C", "T", "Mean"};
	int baseNum=bases.size();
	//
	std::vector<std::vector<uint64_t>> rawBaseQualMerge, rawBaseCountsMerge;
	std::vector<std::vector<uint64_t>> raw3pBaseQualMerge, raw3pBaseCountsMerge;
	rawBaseQualMerge.resize(rawMax, std::vector<uint64_t>(baseNum));
	rawBaseCountsMerge.resize(rawMax, std::vector<uint64_t>(baseNum));
	raw3pBaseQualMerge.resize(150, std::vector<uint64_t>(baseNum));
	raw3pBaseCountsMerge.resize(150, std::vector<uint64_t>(baseNum));

	/////////////////////////////////////////merge data from different threads////////////////////////
	for (int i=0; i<n_threads; i++){
		for (int j=0; j<rawBaseQual[i].size();j++){
			for (int x=0; x<rawBaseQual[i][j].size(); x++){
				rawBaseQualMerge[j][x]+=rawBaseQual[i][j][x];
				rawBaseCountsMerge[j][x]+=rawBaseCounts[i][j][x];
			}
		}
	}
	
	for (int i=0; i<n_threads; i++){
		for (int j=0; j<raw3pBaseQual[i].size();j++){
			for (int x=0; x<raw3pBaseQual[i][j].size(); x++){
				raw3pBaseQualMerge[j][x]+=raw3pBaseQual[i][j][x];
				raw3pBaseCountsMerge[j][x]+=raw3pBaseCounts[i][j][x];
			}
		}
	}
	
    ///////////////////////////////////////////////////whole reads/////////////////////////////////////
    rawReadsQual.x.resize(rawMax);
    rawReadsQual.y.resize(baseNum);
    for (int i = 0; i < baseNum; ++i) {
        rawReadsQual.y[i].bases=bases[i];
        rawReadsQual.y[i].data.resize(rawMax);
    }
    
    rawBasesContents.x.resize(rawMax);
    rawBasesContents.y.resize(baseNum-1);
    for (int i = 0; i < baseNum-1; ++i) {
        rawBasesContents.y[i].bases=bases[i];
        rawBasesContents.y[i].data.resize(rawMax);
    }

	for (int i=0; i < rawMax; i++){
		rawReadsQual.x[i] = i;
        rawBasesContents.x[i] = i;
		rawGCSum += rawBaseCountsMerge[i][1];
        rawGCSum += rawBaseCountsMerge[i][2];
		rawBaseQualSum += rawBaseQualMerge[i][baseNum-1];
		for (int j=0; j<baseNum-1; j++){
			if (rawBaseCountsMerge[i][j] > 0){
				rawReadsQual.y[j].data[i] = rawBaseQualMerge[i][j]/rawBaseCountsMerge[i][j];
			}else{
				rawReadsQual.y[j].data[i] = 0;
			}

			if (rawBaseCountsMerge[i][baseNum-1] > 0){
				rawBasesContents.y[j].data[i]=rawBaseCountsMerge[i][j]/rawBaseCountsMerge[i][baseNum-1];
			}else{
				rawBasesContents.y[j].data[i]=0;
			}
		
		}

		if (rawBaseCountsMerge[i][baseNum-1] >0 ){
			rawReadsQual.y[baseNum-1].data[i]=rawBaseQualMerge[i][baseNum-1]/rawBaseCountsMerge[i][baseNum-1];
		}else{
			rawReadsQual.y[baseNum-1].data[i]=0;
		}
	}
	rawGC = static_cast<float>(rawGCSum*100)/rawBases; // GC content(%)
	//rawGC = round(rawGC * 100) / 100.0;
    rawMeanQual = static_cast<float>(rawBaseQualSum)/rawBases;
	//rawMeanQual = round(rawMeanQual * 100) / 100.0;
	////////////////////////////////////read 5p ////////////////////////////////////////////////////
	raw5pReadsQual.x.resize(150);
    raw5pReadsQual.y.resize(baseNum);
    for (int i = 0; i < baseNum; ++i) {
        raw5pReadsQual.y[i].bases=bases[i];
        raw5pReadsQual.y[i].data.resize(150);
    }
    
    raw5pBasesContents.x.resize(150);
    raw5pBasesContents.y.resize(baseNum-1);
    for (int i = 0; i < baseNum-1; ++i) {
        raw5pBasesContents.y[i].bases=bases[i];
        raw5pBasesContents.y[i].data.resize(150);
    }

	for (int i=0; i < 150; i++){
		raw5pReadsQual.x[i] = i;
        raw5pBasesContents.x[i] = i;
		for (int j=0; j<baseNum-1; j++){
			if (rawBaseCountsMerge[i][j] > 0){
				raw5pReadsQual.y[j].data[i]=rawBaseQualMerge[i][j]/rawBaseCountsMerge[i][j];
			}else{
				raw5pReadsQual.y[j].data[i]=0;
			}

			if (rawBaseCountsMerge[i][baseNum-1] > 0){
				raw5pBasesContents.y[j].data[i]=rawBaseCountsMerge[i][j]/rawBaseCountsMerge[i][baseNum-1];
			}else{
				raw5pBasesContents.y[j].data[i]=0;
			}    
		}

		if (rawBaseCountsMerge[i][baseNum-1] > 0){
			raw5pReadsQual.y[baseNum-1].data[i]=rawBaseQualMerge[i][baseNum-1]/rawBaseCountsMerge[i][baseNum-1];
		}else{
			raw5pReadsQual.y[baseNum-1].data[i]=0;
		}
		
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////read 3p///////////////////////////////////////////////
	
	raw3pReadsQual.x.resize(150);
    raw3pReadsQual.y.resize(baseNum);
    for (int i = 0; i < baseNum; ++i) {
        raw3pReadsQual.y[i].bases=bases[i];
        raw3pReadsQual.y[i].data.resize(150);
    }
    
    raw3pBasesContents.x.resize(150);
    raw3pBasesContents.y.resize(baseNum-1);
    for (int i = 0; i < baseNum-1; ++i) {
        raw3pBasesContents.y[i].bases=bases[i];
        raw3pBasesContents.y[i].data.resize(150);
    }

	for (int i=0; i < 150; i++){
		raw3pReadsQual.x[i] = i;
        raw3pBasesContents.x[i] = i;
		for (int j=0; j<baseNum-1; j++){
			if (raw3pBaseCountsMerge[i][j] > 0){
				raw3pReadsQual.y[j].data[i]=raw3pBaseQualMerge[i][j]/raw3pBaseCountsMerge[i][j];
			}else{
				raw3pReadsQual.y[j].data[i]=0;
			}
			
			if (raw3pBaseCountsMerge[i][baseNum-1]>0){
				raw3pBasesContents.y[j].data[i]=raw3pBaseCountsMerge[i][j]/raw3pBaseCountsMerge[i][baseNum-1];
			}else{
				raw3pBasesContents.y[j].data[i]=0;
			}
		}

		if (raw3pBaseCountsMerge[i][baseNum-1] >0){
			raw3pReadsQual.y[baseNum-1].data[i]=raw3pBaseQualMerge[i][baseNum-1]/raw3pBaseCountsMerge[i][baseNum-1];
		}else{
			raw3pReadsQual.y[baseNum-1].data[i]=0;
		}
		
	}
}

void Get_plot_bar_data(std::vector<int> &rawLens,
						BarPlotData &rawLenDis) {
	int rawMax=rawLens.back();
	std::vector<int> rawLenDisData(int(rawMax/100));
	for (int len : rawLens) {
		int index = int(len/100);
		rawLenDisData[index]++;

	}
	rawLenDis.x.resize(rawLenDisData.size());
	rawLenDis.y.resize(rawLenDisData.size());
	for (int i = 0; i < rawLenDisData.size(); ++i) {
		rawLenDis.x[i] = i*100;
		rawLenDis.y[i] = rawLenDisData[i];
	}	
}


int Get_N50(std::vector<int> &rawLens, int &rawNum, uint64_t &rawBases){
	uint64_t rawN50Bases=0;
	for (int i = rawNum - 1; i >= 0; --i) {
		rawN50Bases += rawLens[i];
		if (rawN50Bases >= rawBases/2) {
			return rawLens[i];
		}
	}
	return 0;
}

std::string limitDecimalPlaces(double num, int places) {
    std::string str = std::to_string(num);

    size_t decimalPos = str.find('.');

    if (decimalPos != std::string::npos && str.size() - decimalPos > (places + 1)) {
        str = str.substr(0, decimalPos + places + 1);
    }

    return str;
}


/// ////////////////////////////////////////////////////////////////
int main (int argc, char *argv[ ]) {
	Para_A24 * P2In = new Para_A24;
	int Inflag=0;
	Inflag=TGSFilter_cmd(argc, argv, P2In);
	if(Inflag==1) {
		delete P2In ;
		return 1 ;
	}

	for (int i=0; i<256;i++) {
		complement[i]='N';
	}
	
	complement['A']='T'; complement['G']='C';  
	complement['C']='G'; complement['T']='A';
	complement['a']='t'; complement['g']='c';  
	complement['c']='g'; complement['t']='a';
	complement['M']='K'; complement['R']='Y';  
	complement['W']='W'; complement['S']='S';  
	complement['Y']='R'; complement['K']='M';
	complement['m']='k'; complement['r']='y';  
	complement['w']='w'; complement['s']='s';  
	complement['y']='r'; complement['k']='m';

	string InPath=(P2In->InFile);
	P2In->Infq = GetFileType(InPath);
	string prefix=GetFilePreifx(InPath);
	string htmlFileName=prefix + ".html";

	string OutPath;
	if (!(P2In->OutFile).empty()) {
		OutPath=(P2In->OutFile);
		htmlFileName = GetFilePreifx(OutPath) + ".html";
		P2In->Outfq = GetFileType(OutPath);
		string ext = GetFileExtension(OutPath);
		if (ext=="gz"){
			P2In->OUTGZ=true;
		}
	}else{
		if (P2In->FastaOut){
			P2In->Outfq=0;
		}else{
			P2In->Outfq = P2In->Infq;
		}
	}

	if ((P2In->Infq)==3 || (P2In->Outfq)==3) {
		cerr<<"Error: The file name suffix should be '.[fastq|fq|fasta|fa][.gz] or .[sam|bam]'"<<endl;
		if ((P2In->Infq)==3) {
			cerr<<"Error: Please check your input file name: "<<(P2In->InFile)<<endl;
		}
		else if((P2In->Outfq)==3) {
			cerr<<"Error: Please check your output file name: "<<(P2In->OutFile)<<endl;
		}else if ((P2In->Outfq)==2){
			cerr<<"Error: Output file only can be fastq or fasta format: "<<(P2In->OutFile)<<endl;
		}
		return 1;
	}else if ((P2In->Infq)==0 && (P2In->Outfq)==1){
		cerr<<"Error: Fasta format input file can't output fastq format file"<<endl;
		return 1;
	}

	//
	std::unordered_map<std::string, int> seqLens;
	std::vector<int> rawLens;
	std::vector<int> cleanLens;
	std::vector<int> downLens;
	uint64_t rawBases;
	uint64_t cleanBases;
	uint64_t downBases;

	std::vector<std::vector<std::string>> tabInfo(29, std::vector<std::string>(3,"0"));
	BarPlotData rawLenDis, cleanLenDis, downLenDis;
	LinePlotData rawReadsQual, raw5pReadsQual, raw3pReadsQual;
	LinePlotData cleanReadsQual, clean5pReadsQual, clean3pReadsQual;
	LinePlotData downReadsQual, down5pReadsQual, down3pReadsQual;
	LinePlotData rawBasesContents, raw5pBasesContents,raw3pBasesContents;
	LinePlotData cleanBasesContents, clean5pBasesContents, clean3pBasesContents;
	LinePlotData downBasesContents, down5pBasesContents, down3pBasesContents;

	string downInput;
	if (P2In->Filter){
		Get_filter_parameter(P2In);
		if (P2In->ONLYAD){
			return 0;
		}

		if (P2In->Downsample){
			string rand=generateRandomString(5);
			if ((P2In->Outfq)==0){
				P2In->TmpOutFile = prefix + ".tmp." + rand + ".fa";
			}else if ((P2In->Outfq)==1){
				P2In->TmpOutFile = prefix + ".tmp."+ rand + ".fq";
			}
		}

		TGSFilterTask task(P2In);
    	task.start();
		downInput=P2In->TmpOutFile;
		seqLens=task.seqLens;
		//
		/////////////////////////////////////raw table//////////////////////////////////////////
		
		rawLens=task.rawLens;
		rawBases=task.rawBases;
		int rawNum=rawLens.size();
		tabInfo[0][0] = std::to_string(rawNum); // reads number befor filtering
		tabInfo[1][0] = std::to_string(rawBases); // bases number befor filtering
		std::sort(rawLens.begin(), rawLens.end());
		int rawMin=rawLens.front();
		int rawMax=rawLens.back();
		tabInfo[3][0] = std::to_string(rawMin); // min
		tabInfo[4][0] = std::to_string(rawMax); // max
		tabInfo[5][0] = std::to_string(int(rawBases/rawNum)); // mean
		tabInfo[6][0] = std::to_string(rawLens[int(rawNum/2)]); // median
		tabInfo[7][0] = std::to_string(Get_N50(rawLens, rawNum, rawBases)); // N50
		/////////////////////////////////////clean table//////////////////////////////////////////
		
		cleanBases=task.cleanBases;
		cleanLens=task.cleanLens;
		int cleanNum=cleanLens.size();
		tabInfo[0][1] = std::to_string(cleanNum); // reads number after filtering
		tabInfo[1][1] = std::to_string(cleanBases); // bases number after filtering
		std::sort(cleanLens.begin(), cleanLens.end());
		int cleanMin=cleanLens.front();
		int cleanMax=cleanLens.back();
		tabInfo[3][1] = std::to_string(cleanMin);
		tabInfo[4][1] = std::to_string(cleanMax);
		tabInfo[5][1] = std::to_string(int(cleanBases/cleanNum));
		tabInfo[6][1] = std::to_string(cleanLens[int(cleanNum/2)]);
		tabInfo[7][1] = std::to_string(Get_N50(cleanLens, cleanNum, cleanBases));
		
		////////////////////////////////////////////////raw reads base and qual plot/////////////////////////////////////////
		
		float rawGC;
		float rawMeanQual;
		Get_plot_line_data(P2In, rawMax, rawGC, rawMeanQual, rawBases,
                    	rawReadsQual, rawBasesContents,
						raw5pReadsQual, raw5pBasesContents,
						raw3pReadsQual, raw3pBasesContents,
                    	task.rawBaseQual, task.rawBaseCounts,
						task.raw3pBaseQual, task.raw3pBaseCounts);
		tabInfo[2][0] = limitDecimalPlaces(round(rawGC * 1000) / 1000.0, 3); // GC content (% ,string)
		tabInfo[8][0] = limitDecimalPlaces(round(rawMeanQual * 1000) / 1000.0, 3); // mean quality
		////////////////////////////////////////////////clean reads base and qual plot/////////////////////////////////////////
		
		float cleanGC;
		float cleanMeanQual;
		Get_plot_line_data(P2In, cleanMax, cleanGC, cleanMeanQual, cleanBases,
                    	cleanReadsQual, cleanBasesContents,
						clean5pReadsQual, clean5pBasesContents,
						clean3pReadsQual, clean3pBasesContents,
                    	task.cleanBaseQual, task.cleanBaseCounts,
						task.clean3pBaseQual, task.clean3pBaseCounts);
		tabInfo[2][1] = limitDecimalPlaces(round(cleanGC * 1000) / 1000.0, 3); // GC content (% ,string)
		tabInfo[8][1] = limitDecimalPlaces(round(cleanMeanQual * 1000) / 1000.0, 3); // mean quality
		/////////////////////////////////////Length distribution/////////////////////////////////////////////
		
		Get_plot_bar_data(rawLens, rawLenDis); // raw Length distribution
		Get_plot_bar_data(cleanLens, cleanLenDis); // clean Length distribution

		/////////////////////////////////raw reads with diffeerent quals//////////////////////
		
		std::vector<uint64_t> rawDiffQualReads;
		Get_diff_quals(task.rawDiffQualReads, rawDiffQualReads);
		tabInfo[9][0]=std::to_string(rawDiffQualReads[0]); // Q5
		tabInfo[10][0]=std::to_string(rawDiffQualReads[1]); //Q7
		tabInfo[11][0]=std::to_string(rawDiffQualReads[2]); //Q10
		tabInfo[12][0]=std::to_string(rawDiffQualReads[3]); //Q13
		tabInfo[13][0]=std::to_string(rawDiffQualReads[4]);	//Q15
		tabInfo[14][0]=std::to_string(rawDiffQualReads[5]); //Q17
		tabInfo[15][0]=std::to_string(rawDiffQualReads[6]); //Q20
		tabInfo[16][0]=std::to_string(rawDiffQualReads[7]); //Q23
		tabInfo[17][0]=std::to_string(rawDiffQualReads[8]); //Q30
		tabInfo[18][0]=std::to_string(rawDiffQualReads[9]); //Q40
		/////////////////////////////////clean reads with diffeerent quals//////////////////////

		std::vector<uint64_t> cleanDiffQualReads;
		Get_diff_quals(task.cleanDiffQualReads, cleanDiffQualReads);
		tabInfo[9][1]=std::to_string(cleanDiffQualReads[0]); // Q5
		tabInfo[10][1]=std::to_string(cleanDiffQualReads[1]); //Q7
		tabInfo[11][1]=std::to_string(cleanDiffQualReads[2]); //Q10
		tabInfo[12][1]=std::to_string(cleanDiffQualReads[3]); //Q13
		tabInfo[13][1]=std::to_string(cleanDiffQualReads[4]); //Q15
		tabInfo[14][1]=std::to_string(cleanDiffQualReads[5]); //Q17
		tabInfo[15][1]=std::to_string(cleanDiffQualReads[6]); //Q20
		tabInfo[16][1]=std::to_string(cleanDiffQualReads[7]); //Q23
		tabInfo[17][1]=std::to_string(cleanDiffQualReads[8]); //Q30
		tabInfo[18][1]=std::to_string(cleanDiffQualReads[9]); //Q40
		
		/////////////////////////////////raw Bases with diffeerent quals//////////////////////

		std::vector<uint64_t> rawDiffQualBases;
		Get_diff_quals(task.rawDiffQualBases, rawDiffQualBases);
		
		tabInfo[19][0]=std::to_string(rawDiffQualBases[0]); // Q5
		tabInfo[20][0]=std::to_string(rawDiffQualBases[1]); //Q7
		tabInfo[21][0]=std::to_string(rawDiffQualBases[2]); //Q10
		tabInfo[22][0]=std::to_string(rawDiffQualBases[3]); //Q13
		tabInfo[23][0]=std::to_string(rawDiffQualBases[4]);	//Q15
		tabInfo[24][0]=std::to_string(rawDiffQualBases[5]); //Q17
		tabInfo[25][0]=std::to_string(rawDiffQualBases[6]); //Q20
		tabInfo[26][0]=std::to_string(rawDiffQualBases[7]); //Q23
		tabInfo[27][0]=std::to_string(rawDiffQualBases[8]); //Q30
		tabInfo[28][0]=std::to_string(rawDiffQualBases[9]); //Q40
		
		/////////////////////////////////clean Bases with diffeerent quals//////////////////////

		std::vector<uint64_t> cleanDiffQualBases;
		Get_diff_quals(task.cleanDiffQualBases, cleanDiffQualBases);
		tabInfo[19][1]=std::to_string(cleanDiffQualBases[0]); // Q5
		tabInfo[20][1]=std::to_string(cleanDiffQualBases[1]); //Q7
		tabInfo[21][1]=std::to_string(cleanDiffQualBases[2]); //Q10
		tabInfo[22][1]=std::to_string(cleanDiffQualBases[3]); //Q13
		tabInfo[23][1]=std::to_string(cleanDiffQualBases[4]); //Q15
		tabInfo[24][1]=std::to_string(cleanDiffQualBases[5]); //Q17
		tabInfo[25][1]=std::to_string(cleanDiffQualBases[6]); //Q20
		tabInfo[26][1]=std::to_string(cleanDiffQualBases[7]); //Q23
		tabInfo[27][1]=std::to_string(cleanDiffQualBases[8]); //Q30
		tabInfo[28][1]=std::to_string(cleanDiffQualBases[9]); //Q40
		///////////////////////////////////////////////////////下面此处数据不要管//////////////////////////////////////////////////////////////////////////////
		uint64_t raw_reads = task.DropInfo[0];
		uint64_t raw_bases = task.DropInfo[1];
		uint64_t clean_reads = task.DropInfo[2];
		uint64_t clean_bases = task.DropInfo[3];
		uint64_t shortDrop_reads = task.DropInfo[4];
		uint64_t shortDrop_bases = task.DropInfo[5];
		uint64_t endDrop_reads = task.DropInfo[6];
		uint64_t endDrop_bases = task.DropInfo[7];
		uint64_t adapterDrop_reads = task.DropInfo[8];
		uint64_t adapterDrop_bases = task.DropInfo[9];
		uint64_t LQDrop_reads = task.DropInfo[10];
		uint64_t LQDrop_bases = task.DropInfo[11];
		uint64_t outDrop_reads = task.DropInfo[12];
		uint64_t outDrop_bases = task.DropInfo[13];
		uint64_t middle_reads = task.DropInfo[14];
		uint64_t middle_bases = task.DropInfo[15];

		cerr << "INFO: "<< raw_reads<<" reads with a total of "<<raw_bases<<" bases were input."<<endl;
		cerr << "INFO: "<< shortDrop_reads <<" reads were discarded with "<<shortDrop_bases<<" bases before filtering."<<endl;
		cerr << "INFO: "<< endDrop_reads <<" reads were trimmed by "<<endDrop_bases<<" bases at the end."<<endl;
		cerr << "INFO: "<< LQDrop_reads <<" reads were discarded with "<<LQDrop_bases<<" bases due to low quality."<<endl;
		cerr << "INFO: "<< middle_reads <<" read was discarded with " <<middle_bases<<" bases due to a middle adapter."<<endl;
		cerr << "INFO: "<< adapterDrop_reads <<" reads were trimmed by " <<adapterDrop_bases<<" bases with an adapter."<<endl;
		cerr << "INFO: "<< outDrop_reads<<" reads were discarded with "<<outDrop_bases<<" bases before output."<<endl;
		cerr << "INFO: "<< clean_reads <<" reads with a total of "<<clean_bases<<" bases were output"<<endl;
		///////////////////////////////////////////////////////上面此处数据不要管//////////////////////////////////////////////////////////////////////////////

	} else {
		downInput=P2In->InFile;
	}

	if (P2In->Downsample){
		DownSampleTask task(P2In, downInput, seqLens);
    	task.start();
		////////////////////////////////////////down table/////////////////////////////////////////
		downBases=task.downBases;
		downLens=task.downLens;
		int downNum=downLens.size();
		tabInfo[0][2] = std::to_string(downNum); // reads number after filtering
		tabInfo[1][2] = std::to_string(downBases); // bases number after filtering
		std::sort(downLens.begin(), downLens.end());
		int downMin=downLens.front();
		int downMax=downLens.back();
		tabInfo[3][2] = std::to_string(downMin);
		tabInfo[4][2] = std::to_string(downMax);
		tabInfo[5][2] = std::to_string(int(downBases/downNum));
		tabInfo[6][2] = std::to_string(downLens[int(downNum/2)]);
		tabInfo[7][2] = std::to_string(Get_N50(downLens, downNum, downBases));
		//
		float downGC;
		float downMeanQual;
		Get_plot_line_data(P2In, downMax, downGC, downMeanQual, downBases,
                    	downReadsQual, downBasesContents,
						down5pReadsQual, down5pBasesContents,
						down3pReadsQual, down3pBasesContents,
                    	task.downBaseQual, task.downBaseCounts,
						task.down3pBaseQual, task.down3pBaseCounts);
		tabInfo[2][2] = limitDecimalPlaces(round(downGC * 1000) / 1000.0, 3); // GC content (% ,string)
		tabInfo[8][2] = limitDecimalPlaces(round(downMeanQual * 1000) / 1000.0, 2); // mean quality
		////////////////////////////////down Length distribution////////////////////////////////////
		Get_plot_bar_data(downLens, downLenDis);

		/////////////////////////////////down reads with diffeerent quals//////////////////////

		std::vector<uint64_t> downDiffQualReads;
		Get_diff_quals(task.downDiffQualReads, downDiffQualReads);
		tabInfo[9][2]=std::to_string(downDiffQualReads[0]); // Q5
		tabInfo[10][2]=std::to_string(downDiffQualReads[1]); //Q7
		tabInfo[11][2]=std::to_string(downDiffQualReads[2]); //Q10
		tabInfo[12][2]=std::to_string(downDiffQualReads[3]); //Q13
		tabInfo[13][2]=std::to_string(downDiffQualReads[4]); //Q15
		tabInfo[14][2]=std::to_string(downDiffQualReads[5]); //Q17
		tabInfo[15][2]=std::to_string(downDiffQualReads[6]); //Q20
		tabInfo[16][2]=std::to_string(downDiffQualReads[7]); //Q23
		tabInfo[17][2]=std::to_string(downDiffQualReads[8]); //Q30
		tabInfo[18][2]=std::to_string(downDiffQualReads[9]); //Q40
		
		/////////////////////////////////down Bases with diffeerent quals//////////////////////

		std::vector<uint64_t> downDiffQualBases;
		Get_diff_quals(task.downDiffQualBases, downDiffQualBases);
		tabInfo[19][2]=std::to_string(downDiffQualBases[0]); // Q5
		tabInfo[20][2]=std::to_string(downDiffQualBases[1]); //Q7
		tabInfo[21][2]=std::to_string(downDiffQualBases[2]); //Q10
		tabInfo[22][2]=std::to_string(downDiffQualBases[3]); //Q13
		tabInfo[23][2]=std::to_string(downDiffQualBases[4]); //Q15
		tabInfo[24][2]=std::to_string(downDiffQualBases[5]); //Q17
		tabInfo[25][2]=std::to_string(downDiffQualBases[6]); //Q20
		tabInfo[26][2]=std::to_string(downDiffQualBases[7]); //Q23
		tabInfo[27][2]=std::to_string(downDiffQualBases[8]); //Q30
		tabInfo[28][2]=std::to_string(downDiffQualBases[9]); //Q30
	}


	///////////////////////////////////////////////////////out test//////////////////////////////////////
	cerr <<"test\tbefore\tafter\tdown"<<endl;
	cerr <<"Total reads\t"<<tabInfo[0][0]<<"\t"<<tabInfo[0][1]<<"\t"<<tabInfo[0][2]<<endl;
	cerr <<"Total bases\t"<<tabInfo[1][0]<<"\t"<<tabInfo[1][1]<<"\t"<<tabInfo[1][2]<<endl;
	cerr <<"GC content(%)\t"<<tabInfo[2][0]<<"\t"<<tabInfo[2][1]<<"\t"<<tabInfo[2][2]<<endl;
	cerr <<"Min Length\t"<<tabInfo[3][0]<<"\t"<<tabInfo[3][1]<<"\t"<<tabInfo[3][2]<<endl;
	cerr <<"Max Length\t"<<tabInfo[4][0]<<"\t"<<tabInfo[4][1]<<"\t"<<tabInfo[4][2]<<endl;
	cerr <<"Mean Length\t"<<tabInfo[5][0]<<"\t"<<tabInfo[5][1]<<"\t"<<tabInfo[5][2]<<endl;
	cerr <<"Median Length\t"<<tabInfo[6][0]<<"\t"<<tabInfo[6][1]<<"\t"<<tabInfo[6][2]<<endl;
	cerr <<"N50 Length\t"<<tabInfo[7][0]<<"\t"<<tabInfo[7][1]<<"\t"<<tabInfo[7][2]<<endl;
	cerr <<"Mean quality\t"<<tabInfo[8][0]<<"\t"<<tabInfo[8][1]<<"\t"<<tabInfo[8][2]<<endl;
	cerr <<"Reads"<<endl;
	cerr <<"Q5\t"<<tabInfo[9][0]<<"\t"<<tabInfo[9][1]<<"\t"<<tabInfo[9][2]<<endl;
	cerr <<"Q7\t"<<tabInfo[10][0]<<"\t"<<tabInfo[10][1]<<"\t"<<tabInfo[10][2]<<endl;
	cerr <<"Q10\t"<<tabInfo[11][0]<<"\t"<<tabInfo[11][1]<<"\t"<<tabInfo[11][2]<<endl;
	cerr <<"Q13\t"<<tabInfo[12][0]<<"\t"<<tabInfo[12][1]<<"\t"<<tabInfo[12][2]<<endl;
	cerr <<"Q15\t"<<tabInfo[13][0]<<"\t"<<tabInfo[13][1]<<"\t"<<tabInfo[13][2]<<endl;
	cerr <<"Q17\t"<<tabInfo[14][0]<<"\t"<<tabInfo[14][1]<<"\t"<<tabInfo[14][2]<<endl;
	cerr <<"Q20\t"<<tabInfo[15][0]<<"\t"<<tabInfo[15][1]<<"\t"<<tabInfo[15][2]<<endl;
	cerr <<"Q23\t"<<tabInfo[16][0]<<"\t"<<tabInfo[16][1]<<"\t"<<tabInfo[16][2]<<endl;
	cerr <<"Q30\t"<<tabInfo[17][0]<<"\t"<<tabInfo[17][1]<<"\t"<<tabInfo[17][2]<<endl;
	cerr <<"Q40\t"<<tabInfo[18][0]<<"\t"<<tabInfo[18][1]<<"\t"<<tabInfo[18][2]<<endl;
	cerr <<"Bases"<<endl;
	cerr <<"Q5\t"<<tabInfo[19][0]<<"\t"<<tabInfo[19][1]<<"\t"<<tabInfo[19][2]<<endl;
	cerr <<"Q7\t"<<tabInfo[20][0]<<"\t"<<tabInfo[20][1]<<"\t"<<tabInfo[20][2]<<endl;
	cerr <<"Q10\t"<<tabInfo[21][0]<<"\t"<<tabInfo[21][1]<<"\t"<<tabInfo[21][2]<<endl;
	cerr <<"Q13\t"<<tabInfo[22][0]<<"\t"<<tabInfo[22][1]<<"\t"<<tabInfo[22][2]<<endl;
	cerr <<"Q15\t"<<tabInfo[23][0]<<"\t"<<tabInfo[23][1]<<"\t"<<tabInfo[23][2]<<endl;
	cerr <<"Q17\t"<<tabInfo[24][0]<<"\t"<<tabInfo[24][1]<<"\t"<<tabInfo[24][2]<<endl;
	cerr <<"Q20\t"<<tabInfo[25][0]<<"\t"<<tabInfo[25][1]<<"\t"<<tabInfo[25][2]<<endl;
	cerr <<"Q23\t"<<tabInfo[26][0]<<"\t"<<tabInfo[26][1]<<"\t"<<tabInfo[26][2]<<endl;
	cerr <<"Q30\t"<<tabInfo[27][0]<<"\t"<<tabInfo[27][1]<<"\t"<<tabInfo[27][2]<<endl;
	cerr <<"Q40\t"<<tabInfo[28][0]<<"\t"<<tabInfo[28][1]<<"\t"<<tabInfo[28][2]<<endl;
	cerr <<"test end"<<endl;
	//////////////////////////////////////////////////////////////////////////////////////////////

	delete P2In ;
	return 0;
}

