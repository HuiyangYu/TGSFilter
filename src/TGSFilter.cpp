#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <thread>
#include <algorithm>
#include <unistd.h>
#include "zlib.h"
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
		"   -l	<int>   min length of read to out [1000]\n"
		"   -q	<int>   min mean base quality [auto]\n"
		"   -5	<int>   drop bases from the front (5') of a read [0]\n"
		"   -3	<int>   drop bases from the tail (3') of a read [0]\n"
		"   -w	<int>   windows szie to cut off low quality region [50]\n"
		"   -t           number of threads [3]\n"
		"   -h           show help [v1.06]\n"
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

int n_thread=3;
int VECMAX =512000; //bp 
int BATCH_SIZE = 24;
int BinWind = 24;

class Para_A24 {
	public:
		int ReadLength;
		int MinLength;
		int AverQ;
		int HeadCrop;
		int TailCrop;
		int Window;
		int Infq;
		int Outfq;
		bool OUTGZ;
		string InFile;
		string OutFile;

		Para_A24() {
			ReadLength=0;
			MinLength=1000;
			AverQ=-1;
			HeadCrop=0;
			TailCrop=0;
			Window=50;
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
			P2In->AverQ=atoi(argv[i]);
		}
		else if (flag == "l") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinLength=atoi(argv[i]);
			if (P2In->MinLength<1){
				P2In->MinLength=1;
			}
		}

		//Trim header and reads options
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
		else if (flag == "w"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->Window=atoi(argv[i]);
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
  	uint32_t len;

	fp = gzopen((P2In->InFile).c_str(), "r");
  	seq = kseq_init(fp);

	int seqNum=0;
	int maxSeq=5000;
	int minQ=50000;
	int maxQ=0;
	uint64_t sumQ = 0;
	uint64_t sumLen = 0;
	int Lengths[maxSeq];
	string qual;
	for (int A=0 ; ((len = kseq_read(seq)) >= 0); A++){
		if (seqNum>=maxSeq){
			break;
		}

		int length=seq->seq.l;
		if (length<10000){
			continue;
		}
		Lengths[seqNum]=length;
		seqNum++;
		qual=seq->qual.s;
		sumLen+=length;
		for (int i=0; i<length; i++) {
			sumQ+=qual[i];
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
	//cout<<"INFO: middle read length: "<<(P2In->ReadLength)<<" bp"<<endl;

	int qType=0;
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
	int meanQ=int(sumQ/sumLen)-qType;
	//cout <<"sumQ: "<<sumQ<<" sumLen: "<<sumLen<<endl;
	//cout<<"INFO: mean base quality: "<<meanQ<<endl;
	if(P2In->AverQ<0){
		if (meanQ>=25){
			P2In->AverQ=20;
		}else if(meanQ>=15){
			P2In->AverQ=10;
		}else {
			P2In->AverQ=7;
		}	
	}
	cout<<"INFO: min mean base quality was set to: "<<(P2In->AverQ)<<endl;
	//cout << "INFO: Qtype "<<qType<<endl;
	P2In->AverQ=(P2In->AverQ)+qType;

	return 0;
}

double CalcAvgQuality(const string &qual) {
	long long int  n = qual.length();
	if (n == 0) {
		return 0;
	}
	long long int sumQ = 0;
	for (long long int i = 0; i < n; i++) {
		sumQ += qual[i];
	}
	return  (sumQ) / n;
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
				string seq = SEQ[i].substr(headCrop, seqLen - totalCrop);
				string qual = QUAL[i].substr(headCrop, qualLen - totalCrop);
				seqLen = seq.length();
				if (seqLen >= (P2In->MinLength)) {
					string out_seq;
					string out_qual;
					string segment_seq;
					string segment_qual;
					double avgQuality;
					string seq_id=ID[i];
					int window=P2In->Window;
					int start_idx = 0;
					int end_idx = window;
					int pass_num=0;
					while (end_idx < seqLen) {
        				segment_seq = seq.substr(start_idx, window);
        				segment_qual = qual.substr(start_idx, window);
						avgQuality = CalcAvgQuality(segment_qual);
						if (avgQuality >= (P2In->AverQ)) {
							out_seq+=segment_seq;
							out_qual+=segment_qual;
						}else{
							if (out_seq.length()>=(P2In->MinLength)){
								pass_num++;
								if (pass_num>=2){
									seq_id=ID[i]+":"+std::to_string(pass_num);
								}
								OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + out_seq + 
										 "\n+\n" + out_qual + "\n";
								out_seq="";
								out_qual="";
							}
						}
						start_idx+=window;
						end_idx=start_idx+window;
					}
					//Processing the final segment
					segment_seq = seq.substr(start_idx);
					segment_qual = qual.substr(start_idx);
					avgQuality = CalcAvgQuality(segment_qual);
					if (avgQuality >= (P2In->AverQ)){
						out_seq+=segment_seq;
						out_qual+=segment_qual;
					}
					pass_num++;
					if (pass_num>=2){
						seq_id=ID[i]+":"+std::to_string(pass_num);
					}
					OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + out_seq + 
								"\n+\n" + out_qual + "\n";
				}
			}
		}
	}else{
		for (int i = Start; i < End; i++) {
			int seqLen = SEQ[i].length();
			int qualLen = QUAL[i].length();
			if ((totalCrop < seqLen) && (seqLen == qualLen)) {
				string seq = SEQ[i].substr(headCrop, seqLen - totalCrop);
				string qual = QUAL[i].substr(headCrop, qualLen - totalCrop);
				seqLen = seq.length();
				if (seqLen >= (P2In->MinLength)) {
					string out_seq;
					string segment_seq;
					string segment_qual;
					double avgQuality;
					string seq_id=ID[i];
					int window=P2In->Window;
					int start_idx = 0;
					int end_idx = window;
					int pass_num=0;
					while (end_idx < seqLen) {
        				segment_seq = seq.substr(start_idx, window);
        				segment_qual = qual.substr(start_idx, window);
						avgQuality = CalcAvgQuality(segment_qual);
						if (avgQuality >= (P2In->AverQ)) {
							out_seq+=segment_seq;
						}else{
							if (out_seq.length()>=(P2In->MinLength)){
								pass_num++;
								if (pass_num>=2){
									seq_id=ID[i]+":"+std::to_string(pass_num);
								}
								OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + out_seq + "\n";
								out_seq="";
							}
						}
						start_idx+=window;
						end_idx=start_idx+window;
					}
					//Processing the final segment
					segment_seq = seq.substr(start_idx);
					segment_qual = qual.substr(start_idx);
					avgQuality = CalcAvgQuality(segment_qual);
					if (avgQuality >= (P2In->AverQ)){
						out_seq+=segment_seq;
					}
					pass_num++;
					if (pass_num>=2){
						seq_id=ID[i]+":"+std::to_string(pass_num);
					}
					OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + out_seq + "\n";
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
			string seq = SEQ[i].substr(headCrop, seqLen - totalCrop);
			seqLen = seq.length();
			if (seqLen >= (P2In->MinLength)) {
				OUT_DATA=OUT_DATA+">"+ID[i]+"\n" + seq + "\n";
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
			int flag=0;
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
				flag++;
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			for (int i = 0; i < flag; i++) {
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
			int flag=0;
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
				flag++;
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			for (int i = 0; i < flag; i++) {
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

