#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <thread>
#include <algorithm>
#include <unistd.h>
#include <zlib.h>
#include "kseq.h"
#include "DeflateOgzstream.cpp"

using namespace std;

KSEQ_INIT(gzFile, gzread)

int  TGSFilter_usage() {
	cout <<""
		"Usage: tgsfilter -1 TGS_reads.fq.gz -o OutFile.fq.gz\n"
		" Options:\n"
		"   -i	<str>   input of fasta/q file\n"
		"   -o	<str>   output of fasta/q file\n"
		"   -q	<int>   min Phred average quality score [10]\n"
		"   -l	<int>   min length of read [100]\n"
		"   -s	<int>   Trim N nucleotides from the start of a read [0]\n"
		"   -e	<int>   Trim N nucleotides from the end of a read [0]\n"
		"   -t           number of threads [1]\n"
		"   -h           show help [v1.05]\n"
		"\n";
	return 1;
}

inline string add_Asuffix (string path) {
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext != "gz") {
		path=path+".gz" ;
	}
	return path ;
}

int n_thread=1;
int VECMAX =512000; //bp 
int BATCH_SIZE = VECMAX;
int BinWind = VECMAX;

class Para_A24 {
	public:
		int minQ;
		int ReadLength;
		int MinLength;
		int HeadCrop;
		int TailCrop;
		int AverQ;
		int Infq;
		int Outfq;
		bool OUTGZ;
		string InFile;
		string OutFile;

		Para_A24() {	
			minQ=10;
			ReadLength=0;
			MinLength=100;
			HeadCrop=0;
			TailCrop=0;
			AverQ=0;
			Infq=2;
			Outfq=2;
			OUTGZ=false;
			InFile="";
			OutFile="";
		}
};

inline void  LogLackArg(string flag) {
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

int TGSFilter_cmd(int argc, char **argv, Para_A24 * P2In) {
	if (argc <= 2) {TGSFilter_usage(); return 1;}

	for(int i = 1; i < argc; i++){
		if(argv[i][0] != '-') {
			cerr << "Error: command option error! please check." << endl;
			return 1;
		}

		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		//input and output options
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

		//Filter short and low quality reads
		else if (flag == "q") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->minQ=atoi(argv[i]);
		}
		else if (flag == "l") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinLength=atoi(argv[i]);
		}

		//Trim header and reads options
		else if (flag == "s") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->HeadCrop=atoi(argv[i]);
		}
		else if (flag == "e"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->TailCrop=atoi(argv[i]);
		}
		//Other options
		else if (flag  ==  "t") {
			if(i + 1 == argc) {LogLackArg(flag) ; return 1;}
			i++;
			n_thread=atoi(argv[i]);
		}
		else if (flag  ==  "u") {
			if(i + 1 == argc) {LogLackArg(flag) ; return 1;}
			i++;
			VECMAX=atoi(argv[i]);
		}
		else if (flag  == "help" || flag  == "h") {
			TGSFilter_usage(); return 1;
		}
		else {
			cerr << "Error: UnKnow argument -"<<flag<<endl;
			return 1;
		}
	}

	// check input and output
	if ((P2In->InFile).empty()) {
		cerr<< "Error: -i lack argument for the must"<<endl;
		return 1;
	}

	if ( (P2In->OutFile).empty()) {
		cerr<< "Error: -o lack argument for the must"<<endl;
		return 1;
	}

	if (access((P2In->InFile).c_str(), 0) != 0) {
		cerr<<"Error: Can't find this file for -i "<<(P2In->InFile)<<endl;
		return 1;
	}
	else {
		return  0;
	}
}

string GetFileExtension(const string& FilePath) {
	size_t dotPos = FilePath.rfind('.');
	if (dotPos == string::npos) {
		return "";
	} 
	else {
		return FilePath.substr(dotPos + 1);
	}
}

int GetFileType(const string& FilePath) {
	int fqfile=2;
	string ext = GetFileExtension(FilePath);
	if (ext == "gz") {
		string FilePathA = FilePath.substr(0, FilePath.rfind('.'));
		ext = GetFileExtension(FilePathA);
	}

	if ((ext == "fq") || (ext == "fastq")) {
		fqfile = 1;
	}
	else if ((ext == "fa") || (ext == "fasta")) {
		fqfile = 0;
	}
	return fqfile;
}

int Get_qType(Para_A24 * P2In){

	gzFile fp;
  	kseq_t *seq;
  	int len;

	fp = gzopen((P2In->InFile).c_str(), "r");
  	seq = kseq_init(fp);

	int seqNum=0;
	int maxSeq=5000;
	int minQ=50000;
	int maxQ=0;
	int Lengths[maxSeq];
	string qual;
	for (int A=0 ; A<(maxSeq) && ((len = kseq_read(seq)) >= 0); A++){
		Lengths[seqNum]=(seq->seq.l);
		seqNum++;
		qual=seq->qual.s;
		for (int i=0; i<(seq->qual.l); i++) {
			if(minQ>qual[i]) {
				minQ=qual[i];
			}
			if(maxQ<qual[i]) {
				maxQ=qual[i];
			}
		}
	}
	kseq_destroy(seq);
  	gzclose(fp);

	sort(Lengths, Lengths + seqNum);
	int middleIndex = seqNum / 2;
	P2In->ReadLength = Lengths[middleIndex];

	int qType=0;
	if (maxQ>0){
		if(minQ >= 33 &&  minQ <= 78  &&  maxQ >= 33 && maxQ <= 78) {
			qType=33;
		}
		else if (minQ >= 64  &&  minQ <= 108  &&  maxQ >= 64 && maxQ <= 108){
			qType=64;
		}
		else if (minQ < 55) {
			qType=33;
		}
		else {
			qType=64;
		}
		P2In->AverQ=(P2In->AverQ)+qType;
	}

	return maxQ;
}

int  CalcAvgQuality(const string& qual) {
	int n = qual.length();
	if (n == 0) {
		return 0;
	}
	long long sumQ = 0;
	for (int i = 0; i < n; i++) {
		sumQ += qual[i];
	}
	return  int (sumQ) / n;
}

void Filter_fastq_reads(Para_A24 * P2In, string &OUT_DATA, int Start, int End, 
						vector<string>& ID, vector <string> &SEQ, vector <string> &QUAL){
	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);

	if ((P2In->Outfq)==1){
		for (int i = Start; i < End; i++) {
			int seqLen = SEQ[i].length();
			int qualLen = QUAL[i].length();
			if ((totalCrop < seqLen) && (seqLen == qualLen)) {
				SEQ[i] = SEQ[i].substr(headCrop, seqLen - totalCrop);
				QUAL[i] = QUAL[i].substr(headCrop, qualLen - totalCrop);
				seqLen = SEQ[i].length();
				if (seqLen >= (P2In->MinLength)) {
					double avgQuality = CalcAvgQuality(QUAL[i]);
					if (avgQuality >= (P2In->AverQ)) {
						OUT_DATA=OUT_DATA+"@"+ID[i]+"\n" + SEQ[i] + "\n+\n" + QUAL[i] + "\n";
					}
				}
			}
		}
	}else{
		for (int i = Start; i < End; i++) {
			int seqLen = SEQ[i].length();
			int qualLen = QUAL[i].length();
			if ((totalCrop < seqLen) && (seqLen == qualLen)) {
				SEQ[i] = SEQ[i].substr(headCrop, seqLen - totalCrop);
				QUAL[i] = QUAL[i].substr(headCrop, qualLen - totalCrop);
				seqLen = SEQ[i].length();
				if (seqLen >= (P2In->MinLength)) {
					double avgQuality = CalcAvgQuality(QUAL[i]);
					if (avgQuality >= (P2In->AverQ)) {
						OUT_DATA=OUT_DATA+">"+ID[i]+"\n" + SEQ[i] + "\n";
					}
				}
			}
		}
	}
}

void Filter_fasta_reads(Para_A24 * P2In,string &OUT_DATA, int Start, int End, 
						vector<string>& ID, vector <string> & SEQ) {
	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);

	for (int i = Start; i < End; i++) {
		int seqLen = SEQ[i].length();
		if (totalCrop < seqLen) {
			SEQ[i] = SEQ[i].substr(headCrop, seqLen - totalCrop);
			seqLen = SEQ[i].length();
			if (seqLen >= (P2In->MinLength)) {
				OUT_DATA=OUT_DATA+">"+ID[i]+"\n" + SEQ[i] + "\n";
			}
		}
	}
}

void Filter_reads(Para_A24 * P2In, string &OUT_DATA, int & Start, int & End, 
				vector<string> &ID, vector <string> &SEQ, vector <string> &QUAL){
	OUT_DATA="";
	if(QUAL[0].length()<=2){
		Filter_fasta_reads(P2In, OUT_DATA, Start, End, ID, SEQ);
	} else {
		Filter_fastq_reads(P2In, OUT_DATA, Start, End, ID, SEQ, QUAL);
	}
}

void Filter_reads_gz(Para_A24 * P2In, int & Start, int & End,
					vector<string> &ID, vector <string> &SEQ, vector <string> &QUAL,
					uint8_t ** ComPresseData, size_t & CompressedSize,
					size_t & OUT_BUFFER_SIZE, int thread_id){
	string OUT_DATA="";
	OUT_DATA.reserve(OUTPUT_BUFFER_SIZE);
	if(QUAL[0].length()<=2){
		Filter_fasta_reads(P2In, OUT_DATA, Start, End, ID, SEQ);

	} else {
		Filter_fastq_reads(P2In, OUT_DATA, Start, End, ID, SEQ, QUAL);
		
	}
	
	if (!(OUT_DATA.empty())) {
		DeflateCompress  GZData;
		size_t inputSize  = OUT_DATA.length();
		GZData.compressData(OUT_DATA.c_str(), inputSize, ComPresseData, 
							CompressedSize,OUT_BUFFER_SIZE,thread_id);
	}
	

}

int Run_seq_filter (Para_A24 * P2In) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  ID;
	vector <string>  SEQ;
	vector <string>  QUAL;
	
	ID.resize(BATCH_SIZE+2);
	SEQ.resize(BATCH_SIZE+2);
	QUAL.resize(BATCH_SIZE+2);

	vector<thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	int len;
	gzFile fp;
	kseq_t *seq;

	fp = gzopen((P2In->InFile).c_str(), "r");
	seq = kseq_init(fp);

	string outname=P2In->OutFile;

	int seq_num=0;

	if (P2In->OUTGZ){
		DeflateOgzstream OUTHGZ(outname.c_str());
		uint8_t ** ComPresseData = new uint8_t*[n_thread];
		size_t * CompressedSize =new size_t [n_thread];
		size_t * OUT_BUFFER_SIZE =new size_t [n_thread];
		int * ArryThread=new int [n_thread];
		for (int i = 0; i < n_thread; i++) {
			ComPresseData[i] = new uint8_t[OUTPUT_BUFFER_SIZE];
			CompressedSize[i]=OUTPUT_BUFFER_SIZE;
			OUT_BUFFER_SIZE[i]=OUTPUT_BUFFER_SIZE;
			ArryThread[i]=i;
		}

		while((len = kseq_read(seq)) >= 0 ) {
			ID[seq_num]=seq->name.s;
			SEQ[seq_num]=seq->seq.s;
			QUAL[seq_num]=seq->qual.s;
			
			seq_num++;
			if (seq_num==BATCH_SIZE) {
				for (int i = 0; i < n_thread; i++) {
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_num){
						End[i]=seq_num;
					}
					if (Start[i]>=End[i]) {
						continue;
					}
					threads.push_back(thread(Filter_reads_gz,P2In,ref(Start[i]),ref(End[i]), 
							ref(ID), ref(SEQ),ref(QUAL),ComPresseData, ref(CompressedSize[i]),
							ref(OUT_BUFFER_SIZE[i]), ref(ArryThread[i])));
				}

				for (auto& thread : threads){
					thread.join();
				}
				threads.clear();

				for (int i = 0; i < n_thread; i++) {
					if (CompressedSize[i]>0) {
						OUTHGZ.writeGZIO(ComPresseData[i], CompressedSize[i] );
					}
				}

				seq_num=0;
			}
		}

		if (seq_num!=0) {
			for (int i = 0; i < n_thread; i++) {
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>seq_num){
					End[i]=seq_num;
				}
				if (Start[i]>=End[i]) {
					continue;
				}
				threads.push_back(thread(Filter_reads_gz,P2In,ref(Start[i]),ref(End[i]), 
						ref(ID), ref(SEQ),ref(QUAL),ComPresseData, ref(CompressedSize[i]),
						ref(OUT_BUFFER_SIZE[i]), ref(ArryThread[i])));
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			for (int i = 0; i < n_thread; i++) {
				if (CompressedSize[i]>0) {
					OUTHGZ.writeGZIO(ComPresseData[i], CompressedSize[i] );
				}
			}

			seq_num=0;
		}

		for (int i = 0; i < n_thread; i++) {
			delete[]  ComPresseData[i] ;
		}
		
		delete [] CompressedSize ;
		delete [] OUT_BUFFER_SIZE;
		delete [] ComPresseData ;
		delete [] ArryThread ;

	} else {
		ofstream OUTH;
		OUTH.open(outname.c_str());
		vector <string> OUT_DATA;
		OUT_DATA.resize(n_thread);

		while((len = kseq_read(seq)) >= 0 ) {
			ID[seq_num]=seq->name.s;
			SEQ[seq_num]=seq->seq.s;
			QUAL[seq_num]=seq->qual.s;
			
			seq_num++;
			if (seq_num==BATCH_SIZE) {
				for (int i = 0; i < n_thread; i++) {
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_num){
						End[i]=seq_num;
					}
					if (Start[i]>=End[i]) {
						continue;
					}
					threads.push_back(thread(Filter_reads,P2In, ref(OUT_DATA[i]), 
					ref(Start[i]),ref(End[i]),ref(ID),ref(SEQ),ref(QUAL)));
				}

				for (auto& thread : threads){
					thread.join();
				}
				threads.clear();

				for (int i = 0; i < n_thread; i++) {
					if (!(OUT_DATA[i].empty())) {
						OUTH << OUT_DATA[i];
					}
				}
				seq_num=0;
			}
		}

		if (seq_num!=0) {
			for (int i = 0; i < n_thread; i++) {
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>seq_num){
					End[i]=seq_num;
				}
				if (Start[i]>=End[i]) {
					continue;
				}
				threads.push_back(thread(Filter_reads,P2In, ref(OUT_DATA[i]), 
				ref(Start[i]),ref(End[i]),ref(ID),ref(SEQ),ref(QUAL)));
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			for (int i = 0; i < n_thread; i++) {
				if (!(OUT_DATA[i].empty())) {
					OUTH << OUT_DATA[i];
				}
			}
			seq_num=0;
		}
		OUTH.close();
	}

	kseq_destroy(seq);
	gzclose(fp);

	delete [] Start;
	delete [] End;

	return 0;
}

int Run_TGSFilter(Para_A24 * P2In) {
	string InPath=(P2In->InFile);
	string OutPath=(P2In->OutFile);	
	P2In->Infq = GetFileType(InPath);
	P2In->Outfq = GetFileType(OutPath);

	if ((P2In->Infq)==2 || (P2In->Outfq)==2) {
		cerr<<"Error: The file name suffix should be '.[fastq|fq|fasta|fa][.gz]'"<<endl;
		if ((P2In->Infq)==2) {
			cerr<<"Error: Please check your input file name: "<<(P2In->InFile)<<endl;
		}
		else if((P2In->Outfq)==2) {
			cerr<<"Error: Please check your output file name: "<<(P2In->OutFile)<<endl;
		}
		return 1;
	}else if ((P2In->Infq)==0 && (P2In->Outfq)==1){
		cerr<<"Fasta format input file can't output fastq format file"<<endl;
		return 1;
	}

	string ext = GetFileExtension(OutPath);
	if (ext=="gz"){
		P2In->OUTGZ=true;
	}

	Get_qType(P2In);
	int read_length=(P2In->ReadLength);
	cout<<"INFO: middle read length is "<<read_length<<" bp"<<endl;
	
	BinWind=int(VECMAX/read_length);
	if (BinWind<1){
		BinWind=1;
	}
	BATCH_SIZE = BinWind*n_thread;

	Run_seq_filter(P2In);

	return 0;
}

//////////////////main///////////////////
int main (int argc, char *argv[ ]) {
	Para_A24 * P2In = new Para_A24;
	int Inflag=0;
	Inflag=TGSFilter_cmd(argc, argv, P2In);
	if(Inflag==1) {
		delete P2In ;
		return 1 ;
	}

	Run_TGSFilter(P2In);

	delete P2In ;
	return 0;
}

