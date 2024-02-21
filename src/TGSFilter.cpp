#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <thread>
#include <algorithm>
#include <unistd.h>
#include "kseq.h"
#include "gzstream.C"
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
		"   -t          number of threads [1]\n"
		"   -h          show help [v1.04]\n"
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
int VECMAX =1024;
int BATCH_SIZE = VECMAX;
int BinWind = VECMAX;

class Para_A24 {
	public:
		int minQ;
		int MinLength;
		int HeadCrop;
		int TailCrop;
		int AverQ;
		string InFile;
		string OutFile;

		Para_A24()
		{	
			minQ=10;
			MinLength=100;
			HeadCrop=0;
			TailCrop=0;
			AverQ=0;
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
			if(i + 1 == argc) { LogLackArg(flag) ; return 1;}
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

pair<int, int> GetFileType(const string& FilePath) {
	int fqfile=2;
	int gzfile=2;
	string ext = GetFileExtension(FilePath);

	if (ext == "gz") {
		gzfile = 1;
		string FilePathA = FilePath.substr(0, FilePath.rfind('.'));
		ext = GetFileExtension(FilePathA);
	}
	else {
		gzfile= 0;
	}

	if ((ext == "fq") || (ext == "fastq")) {
		fqfile = 1;
	}
	else if ((ext == "fa") || (ext == "fasta")) {
		fqfile = 0;
	}
	return make_pair(fqfile, gzfile);
}

int GetReadLen(Para_A24 * P2In, int &Ingz){
	igzstream INH ((P2In->InFile).c_str(),ifstream::in);
	INH.rdbuf()->pubsetbuf(nullptr, 1024*8);
	if (INH.fail()) {
		cerr << "Error: Can't open input file: " << (P2In->InFile) << endl;
		return  -1;
	}

	int seqNum=0;
	int maxSeq=5000;
	int minQ=50000;
	int maxQ=0;
	int Lengths[maxSeq];

	string id, seq;

	for (int A=1; A<maxSeq && (!INH.eof()); A++) {
		getline(INH, id);
		getline(INH, seq);

		if (id.empty()) {continue;}
		seqNum++;
		int seqlLen = seq.length();
		Lengths[seqNum]=seqlLen;
	}

	INH.close();

	sort(Lengths, Lengths + maxSeq);
	int middleIndex = maxSeq / 2;
	int middleValue = Lengths[middleIndex];

	return middleValue;
}

pair<int, int> GetQtype(Para_A24 * P2In, int &Ingz) {

	igzstream INH ((P2In->InFile).c_str(),ifstream::in);
	INH.rdbuf()->pubsetbuf(nullptr, 1024*8);
	if (INH.fail()) {
		cerr << "Error: Can't open input file: " << (P2In->InFile) << endl;
		return  make_pair(-1, -1);;
	}

	int seqNum=0;
	int maxSeq=1000;
	int minQ=50000;
	int maxQ=0;
	int Lengths[maxSeq];

	string id, seq, plus, qual;

	for (int A=1; A<maxSeq && (!INH.eof()); A++)
	{
		getline(INH, id);
		getline(INH, seq);
		getline(INH, plus);
		getline(INH, qual);

		if (id.empty()) {continue;}

		seqNum++;

		int qualLen = qual.length();
		Lengths[seqNum]=qualLen;

		for (int i=0; i<qualLen; i++) {
			if(minQ>qual[i]) {
				minQ=qual[i];
			}
			if(maxQ<qual[i]) {
				maxQ=qual[i];
			}
		}
	}

	INH.close();

	sort(Lengths, Lengths + maxSeq);
	int middleIndex = maxSeq / 2;
	int middleValue = Lengths[middleIndex];

	int qType;
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

	return make_pair(qType, middleValue);
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



void FilterFQFQGZ(Para_A24 * para_, bool & pass_, uint8_t * ComPresseData, size_t & CompressedSize,size_t &  OUT_BUFFER_SIZE , int & start_, int & end_, vector <string> & id_, vector <string> & seq_, vector <string> & qual_){
	int headCrop = para_->HeadCrop;
	int totalCrop = (para_->HeadCrop) + (para_->TailCrop);
	string  Data=""; Data.reserve(OUTPUT_BUFFER_SIZE);

	for (int i = start_; i < end_; i++) {
		int seqLen = seq_[i].length();
		int qualLen = qual_[i].length();
		if ((totalCrop >= seqLen) || (seqLen != qualLen)) {
		}
		else {
			seq_[i] = seq_[i].substr(headCrop, seqLen - totalCrop);
			qual_[i] = qual_[i].substr(headCrop, qualLen - totalCrop);
			seqLen = seq_[i].length();
			if (seqLen < (para_->MinLength)) {
			}
			else {
				double avgQuality = CalcAvgQuality(qual_[i]);
				if (avgQuality < (para_->AverQ)) {
				}
				else {
					Data=Data+id_[i]+"\n" + seq_[i] + "\n+\n" + qual_[i] + "\n";
				}
			}
		}
	}

	if (Data.empty()) {
		pass_ = false;
	}
	else {
		pass_ = true;
		DeflateCompress  GZData ;
		size_t inputSize  = Data.length();
		GZData.compressData(Data.c_str(), inputSize, ComPresseData , CompressedSize,OUT_BUFFER_SIZE);
	}
}

void FilterFQFAGZ(Para_A24  * para_, bool  &  pass_,  uint8_t * ComPresseData , size_t & CompressedSize ,size_t &  OUT_BUFFER_SIZE  , int & start_, int & end_, vector <string> & id_, vector <string> & seq_, vector <string> & qual_)
{
	int headCrop = para_->HeadCrop;
	int totalCrop = (para_->HeadCrop) + (para_->TailCrop);
	string  Data=""; Data.reserve(OUTPUT_BUFFER_SIZE);
	for (int i = start_; i < end_; i++) 
	{
		int seqLen = seq_[i].length();
		int qualLen = qual_[i].length();
		if ((totalCrop >= seqLen) || (seqLen != qualLen)) {
		}
		else {
			seq_[i] = seq_[i].substr(headCrop, seqLen - totalCrop);
			qual_[i] = qual_[i].substr(headCrop, qualLen - totalCrop);
			seqLen = seq_[i].length();
			if (seqLen < (para_->MinLength)) {
			}
			else {
				double avgQuality = CalcAvgQuality(qual_[i]);
				if (avgQuality < (para_->AverQ)) {
				}
				else{
					id_[i][0]='>';
					Data=Data+id_[i]+"\n" + seq_[i] +"\n";
				}
			}
		}
	}

	if (Data.empty()) {
		pass_ = false;
	}
	else {
		pass_ = true;
		DeflateCompress  GZData ;
		size_t inputSize  = Data.length();
		GZData.compressData(Data.c_str(), inputSize, ComPresseData , CompressedSize,OUT_BUFFER_SIZE);
	}
}

void FilterFAFAGZ(Para_A24  * para_, bool  &  pass_,  uint8_t * ComPresseData , size_t & CompressedSize,size_t &  OUT_BUFFER_SIZE , int & start_, int & end_, vector <string> & id_, vector <string> & seq_){
	int headCrop = para_->HeadCrop;
	int totalCrop = (para_->HeadCrop) + (para_->TailCrop);
	string  Data=""; Data.reserve(OUTPUT_BUFFER_SIZE);
	for (int i = start_; i < end_; i++) {
		int seqLen = seq_[i].length();
		if ((totalCrop >= seqLen)) {
		}
		else {
			seq_[i] = seq_[i].substr(headCrop, seqLen - totalCrop);
			seqLen = seq_[i].length();
			if (seqLen < (para_->MinLength)) {
			}
			else {
				id_[i][0]='>';
				Data=Data+id_[i]+"\n" + seq_[i] +"\n";
			}
		}
	}

	if (Data.empty()) {
		pass_ = false;
	}
	else {
		pass_ = true;
		DeflateCompress GZData ;
		size_t inputSize = Data.length();
		GZData.compressData(Data.c_str(), inputSize, ComPresseData , CompressedSize,OUT_BUFFER_SIZE);
	}
}

void FilterFQ(Para_A24  * para_, bool  *  pass_, int & start_, int & end_, vector <string> & seq_, vector <string> & qual_) {
	int headCrop = para_->HeadCrop;
	int totalCrop = (para_->HeadCrop) + (para_->TailCrop);
	for (int i = start_; i < end_; i++) {
		int seqLen = seq_[i].length();
		int qualLen = qual_[i].length();
		if ((totalCrop >= seqLen) || (seqLen != qualLen)) {
			pass_[i] = false;
		}
		else {
			seq_[i] = seq_[i].substr(headCrop, seqLen - totalCrop);
			qual_[i] = qual_[i].substr(headCrop, qualLen - totalCrop);
			seqLen = seq_[i].length();
			if (seqLen < (para_->MinLength)) {
				pass_[i] = false;
			}
			else {
				double avgQuality = CalcAvgQuality(qual_[i]);
				if (avgQuality < (para_->AverQ)) {
					pass_[i] = false;
				} else {
					pass_[i] = true;
				}
			}
		}
	}
}

void FilterFA(Para_A24* para_, bool * pass_, int & start_, int & end_, vector <string> & seq_ ) {
	int headCrop = para_->HeadCrop;
	int totalCrop = (para_->HeadCrop) + (para_->TailCrop);
	int minLength = para_->MinLength;
	for (int i = start_; i < end_; i++) {
		int seqLen = seq_[i].length();
		if (totalCrop >= seqLen) {
			pass_[i] = false;
		}
		else {
			seq_[i] = seq_[i].substr(headCrop, seqLen - totalCrop);
			seqLen = seq_[i].length();
			if (seqLen < minLength) {
				pass_[i] = false;
			}
			else {
				pass_[i] = true;
			}
		}
	}
}

void compressFasta(vector<string>& ID, vector<string>& SEQ, bool * PASS,DeflateOgzstream & OUTH, int end){
	string  input="";
	for (int j = 0; j < end; j++) {
		if (PASS[j]) {
			input=input+ ID[j] + "\n" + SEQ[j] + "\n";
		}
	}
	int length= input.length();
	if (length>1) {
		OUTH.writeGZIO(input);
	}
}

void compressFastq(vector<string>& ID, vector<string>& SEQ, vector<string>& QUAL, bool * PASS, DeflateOgzstream & OUTH, int end){
	string input="";
	for (int j = 0 ; j < end; j++) {
		if (PASS[j]) {
			input=input+ ID[j] + "\n" + SEQ[j] + "\n+\n" + QUAL[j] + "\n";
		}
	}
	int length= input.length();
	if (length>1) {
		OUTH.writeGZIO(input);
	}
}

int Run_fasta(Para_A24 * P2In, int &Ingz, int &Outgz) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  ID;
	vector <string>  SEQ;

	ID.resize(BATCH_SIZE+2);
	SEQ.resize(BATCH_SIZE+2);

	int A=0;

	vector<thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool *PASS =new bool [BATCH_SIZE];

	gzFile FA;
	kseq_t *seqFA;
	FA = gzopen((P2In->InFile).c_str(), "r");
	seqFA = kseq_init(FA);

	string id, seq, plus, qual;
	long long  seq_number=0;

	if (Outgz == 1) {
		P2In->OutFile=add_Asuffix(P2In->OutFile);
		DeflateOgzstream OUTH(P2In->OutFile);
		int AA;
		uint8_t ** ComPresseData = new uint8_t*[n_thread];
		size_t * CompressedSize =new size_t [n_thread];
		size_t * OUT_BUFFER_SIZE =new size_t [n_thread];

		for (int i = 0; i < n_thread; i++) {
			ComPresseData[i] = new uint8_t[OUTPUT_BUFFER_SIZE];
			CompressedSize[i]=OUTPUT_BUFFER_SIZE;
			OUT_BUFFER_SIZE[i]=OUTPUT_BUFFER_SIZE;
		}

		while((AA = kseq_read(seqFA)) >= 0 ) {
			id=(seqFA->name.s);
			seq=(seqFA->seq.s);

			ID[seq_number]=id;
			SEQ[seq_number]=seq;
			seq_number++;

			if (seq_number==BATCH_SIZE) {
				for (int i = 0; i < n_thread; i++) {
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_number) {
						End[i]=seq_number;
					}

					if (Start[i]>=End[i]) {
						continue;
					}

					threads.push_back(thread(FilterFAFAGZ, P2In, ref(PASS[i]),ComPresseData[i], ref(CompressedSize[i]),ref(OUT_BUFFER_SIZE[i]), ref(Start[i]),ref(End[i]), ref(ID), ref(SEQ)));

				}

				for (auto& thread : threads) {
					thread.join();
				}
				threads.clear();
				for (int i = 0; i < n_thread; i++) {
					if (PASS[i]) {
						OUTH.writeGZIO ( ComPresseData[i] , CompressedSize[i] );
					}
				}

				seq_number=0;
			}
		}

		if (seq_number!=0) {
			for (int i = 0; i < n_thread; i++){
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>seq_number){
					End[i]=seq_number;
				}
				if (Start[i]>=End[i]) {
					continue;
				}
				threads.push_back(thread(FilterFAFAGZ, P2In, ref(PASS[i]),ComPresseData[i], ref(CompressedSize[i]),ref(OUT_BUFFER_SIZE[i]), ref(Start[i]),ref(End[i]), ref(ID), ref(SEQ)));
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			for (int i = 0; i < n_thread; i++) {
				if (PASS[i]) {
					OUTH.writeGZIO (ComPresseData[i], CompressedSize[i]);
				}
			}
			seq_number=0;
		}


		for (int i = 0; i < n_thread; i++) {
			delete[]  ComPresseData[i] ;
		}

		delete [] ComPresseData ;
		delete [] CompressedSize ;
	} 
	else {
		ofstream OUTH;
		OUTH.open((P2In->OutFile));

		int AA;
		while((AA = kseq_read(seqFA)) >= 0) {
			id=(seqFA->name.s);
			seq=(seqFA->seq.s);

			ID[seq_number]=id;
			SEQ[seq_number]=seq;

			seq_number++;

			if (seq_number==BATCH_SIZE) {
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_number) {
						End[i]=seq_number;
					}

					if (Start[i]>=End[i]) {
						continue;
					}

					threads.push_back(thread(FilterFA,P2In,PASS,ref(Start[i]),
											 ref(End[i]),ref(SEQ)));
				}

				for (auto& thread : threads) {
					thread.join();
				}
				threads.clear();

				for (int j = 0; j < seq_number; j++) {
					if (PASS[j]) {
						ID[j][0]='>';
						OUTH << ID[j] << "\n" << SEQ[j] << "\n";
					}
				}

				seq_number=0;
			}
		}

		if (seq_number!=0) {
			for (int i = 0; i < n_thread; i++){
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>seq_number){
					End[i]=seq_number;
				}
				if (Start[i]>=End[i]) {
					continue;
				}
				threads.push_back(thread(FilterFA,P2In,PASS,ref(Start[i]),
										 ref(End[i]),ref(SEQ)));
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			for (int j = 0; j < seq_number; j++) {
				if (PASS[j]) {
					ID[j][0]='>';
					OUTH << ID[j] << "\n" << SEQ[j] << "\n";
				}
			}

			seq_number=0;
			OUTH.close();

		}
	}

	kseq_destroy(seqFA);
	gzclose(FA);

	delete [] Start;
	delete [] End;
	delete [] PASS;

	return 0;
}

int Run_fastq(Para_A24 * P2In, int &Ingz, int &Outgz, int &Outfq) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  ID;
	vector <string>  SEQ;
	vector <string>  QUAL;

	ID.resize(BATCH_SIZE+2);
	SEQ.resize(BATCH_SIZE+2);
	QUAL.resize(BATCH_SIZE+2);

	int A=0;

	vector<thread> threads;


	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool *PASS =new bool [BATCH_SIZE];

	igzstream INH ((P2In->InFile).c_str(),ifstream::in);
	INH.rdbuf()->pubsetbuf(nullptr, 1024*8);

	if (INH.fail()) {
		cerr << "Error: Can't open input file: " << (P2In->InFile) << endl;
		return  -1;
	}

	string id, seq, plus, qual;
	long long  seq_number=0;

	//open out file
	if (Outgz == 1) {
		P2In->OutFile=add_Asuffix(P2In->OutFile);
		DeflateOgzstream OUTH(P2In->OutFile);

		uint8_t ** ComPresseData = new uint8_t*[n_thread];
		size_t * CompressedSize =new size_t [n_thread];
		size_t * OUT_BUFFER_SIZE =new size_t [n_thread];
		for (int i = 0; i < n_thread; i++) {
			ComPresseData[i] = new uint8_t[OUTPUT_BUFFER_SIZE];
			CompressedSize[i]=OUTPUT_BUFFER_SIZE;
			OUT_BUFFER_SIZE[i]=OUTPUT_BUFFER_SIZE;
		}

		if (Outfq==1) {
			while(!INH.eof()) {
				getline(INH, id);
				getline(INH, seq);
				getline(INH, plus);
				getline(INH, qual);

				if (id.empty()) { 
					continue;
				}

				ID[seq_number]=id;
				SEQ[seq_number]=seq;
				QUAL[seq_number]=qual;

				seq_number++;

				if (seq_number==BATCH_SIZE) {
					for (int i = 0; i < n_thread; i++) {
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						if (End[i]>seq_number) {
							End[i]=seq_number;
						}

						if (Start[i]>=End[i]) {
							continue;
						}
						threads.push_back(std::thread(FilterFQFQGZ,P2In,ref(PASS[i]),ComPresseData[i],ref(CompressedSize[i]), ref(OUT_BUFFER_SIZE[i]),ref(Start[i]),ref(End[i]),ref(ID),ref(SEQ),ref(QUAL)));
					}
					for (auto& thread : threads) {
						thread.join();
					}
					threads.clear();

					for (int i = 0; i < n_thread; i++) {
						if (PASS[i]) {
							OUTH.writeGZIO ( ComPresseData[i] , CompressedSize[i] );
						}
					}
					seq_number=0;
				}
			}

			if (seq_number!=0) {
				for (int i = 0; i < n_thread; i++) {
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_number){
						End[i]=seq_number;
					}
					if (Start[i]>=End[i]) {
						continue;
					}
					threads.push_back(std::thread(FilterFQFQGZ,P2In,ref(PASS[i]),ComPresseData[i],ref(CompressedSize[i]), ref(OUT_BUFFER_SIZE[i]),ref(Start[i]),ref(End[i]),ref(ID),ref(SEQ),ref(QUAL)));
				}

				for (auto& thread : threads){
					thread.join();
				}
				threads.clear();

				//outfile

				for (int i = 0; i < n_thread; i++) {
					if (PASS[i]) {
						OUTH.writeGZIO ( ComPresseData[i] , CompressedSize[i] );
					}
				}

				seq_number=0;
			}
		}
		else {
			while(!INH.eof()) {
				getline(INH, id);
				getline(INH, seq);
				getline(INH, plus);
				getline(INH, qual);
				if (id.empty()) { 
					continue;
				}

				ID[seq_number]=id;
				SEQ[seq_number]=seq;
				QUAL[seq_number]=qual;

				seq_number++;

				if (seq_number==BATCH_SIZE) {
					for (int i = 0; i < n_thread; i++) {
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						if (End[i]>seq_number) {
							End[i]=seq_number;
						}

						if (Start[i]>=End[i]) {
							continue;
						}
						threads.push_back(thread(FilterFQFAGZ,P2In,ref(PASS[i]),
						ComPresseData[i],ref(CompressedSize[i]),ref(OUT_BUFFER_SIZE[i]),ref(Start[i]),
						ref(End[i]),ref(ID),ref(SEQ),ref(QUAL)));
					}
					for (auto& thread : threads) {
						thread.join();
					}
					threads.clear();
					for (int i = 0; i < n_thread; i++) {
						if (PASS[i]) {
							OUTH.writeGZIO ( ComPresseData[i] , CompressedSize[i] );
						}
					}
					seq_number=0;
				}
			}

			if (seq_number!=0) {
				for (int i = 0; i < n_thread; i++) {
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_number){
						End[i]=seq_number;
					}
					if (Start[i]>=End[i]) {
						continue;
					}
					threads.push_back(thread(FilterFQFAGZ,P2In,ref(PASS[i]),
					ComPresseData[i],ref(CompressedSize[i]),ref(OUT_BUFFER_SIZE[i]),ref(Start[i]),
					ref(End[i]),ref(ID),ref(SEQ),ref(QUAL)));
				}

				for (auto& thread : threads){
					thread.join();
				}
				threads.clear();

				for (int i = 0; i < n_thread; i++) {
					if (PASS[i]) {
						OUTH.writeGZIO ( ComPresseData[i] , CompressedSize[i] );
					}
				}

				seq_number=0;
			}
		}

		for (int i = 0; i < n_thread; i++) {
			delete[]  ComPresseData[i] ;
		}
		
		delete [] CompressedSize ;
		delete [] OUT_BUFFER_SIZE;
		delete [] ComPresseData ;
	} 
	else {	
		ofstream OUTH;
		OUTH.open((P2In->OutFile));
		while(!INH.eof()) {	
			getline(INH, id);
			getline(INH, seq);
			getline(INH, plus);
			getline(INH, qual);
			if (id.empty()) { 
				continue;
			}

			ID[seq_number]=id;
			SEQ[seq_number]=seq;
			QUAL[seq_number]=qual;

			seq_number++;

			if (seq_number==BATCH_SIZE) {
				for (int i = 0; i < n_thread; i++) {
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>seq_number) {
						End[i]=seq_number;
					}

					if (Start[i]>=End[i]) {
						continue;
					}
					threads.push_back(std::thread(FilterFQ,P2In,PASS,ref(Start[i]),
					ref(End[i]),ref(SEQ),ref(QUAL)));
				}

				for (auto& thread : threads) {
					thread.join();
				}
				threads.clear();

				if (Outfq==1) {
					for (int j = 0; j < seq_number; j++) {
						if (PASS[j]) {
							OUTH << ID[j] << "\n" << SEQ[j] << "\n+\n" << QUAL[j] << "\n";
						}
					}
				}
				else if(Outfq==0) {
					for (int j = 0; j < seq_number; j++) {
						if (PASS[j]) {
							ID[j][0]='>';
							OUTH << ID[j] << "\n" << SEQ[j] << "\n";
						}
					}
				}

				seq_number=0;
			}
		}

		if (seq_number!=0) {
			for (int i = 0; i < n_thread; i++){
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>seq_number){
					End[i]=seq_number;
				}
				if (Start[i]>=End[i]) {
					continue;
				}
				threads.push_back(thread(FilterFQ,P2In,PASS,
				ref(Start[i]),ref(End[i]),ref(SEQ),ref(QUAL)));
			}

			for (auto& thread : threads){
				thread.join();
			}
			threads.clear();

			if (Outfq==1) {
				for (int j = 0; j < seq_number; j++) {
					if (PASS[j]) {
						OUTH << ID[j] << "\n" << SEQ[j] << "\n+\n" << QUAL[j] << "\n";
					}
				}
			}
			else if(Outfq==0) {
				for (int j = 0; j < seq_number; j++) {
					if (PASS[j]) {
						ID[j][0]='>';
						OUTH << ID[j] << "\n" << SEQ[j] << "\n";
					}
				}
			}

			seq_number=0;
			OUTH.close();

		}
	}

	INH.close();
	delete [] Start;
	delete [] End;
	delete [] PASS;

	return 0;
}

int Run_TGSFilter(Para_A24 * P2In)
{	
	string InPath=(P2In->InFile);
	pair<int, int> Intype = GetFileType(InPath);
	int Infq = Intype.first;
	int Ingz = Intype.second;

	string OutPath=(P2In->OutFile);
	pair<int, int> Outtype = GetFileType(OutPath);
	int Outfq = Outtype.first;
	int Outgz = Outtype.second;

	if (Infq==2 || Outfq==2) {
		cerr<<"Error: The file name suffix should be '.[fastq|fq|fasta|fa][.gz]'"<<endl;
		if (Infq==2) {
			cerr<<"Error: Please check your input file name: "<<(P2In->InFile)<<endl;
		}
		else if(Outfq==2) {
			cerr<<"Error: Please check your output file name: "<<(P2In->OutFile)<<endl;
		}
		return 1;
	}

	int readLen;
	if (Infq==1) {
		pair<int, int> fqInfo=GetQtype(P2In, Ingz);
		int Qtype=fqInfo.first; // Phred 33 or 64
		P2In->AverQ=(P2In->minQ)+Qtype;
		readLen=fqInfo.second;
	}
	else if (Infq==0) {
		readLen=GetReadLen(P2In, Ingz);
	}

	if (VECMAX == 1024) {
		if (readLen>=1000) {
			VECMAX=8;
		} 
		else
		{
			VECMAX=1024;
		}
	}

	BinWind = int(VECMAX/n_thread)+1;
	if (BinWind<4) {BinWind=4;}
	BATCH_SIZE = BinWind*n_thread;

	if (Infq==1) {
		Run_fastq(P2In, Ingz, Outgz, Outfq);
		cerr<<"done Run FQ"<<endl;
	}
	else if (Infq==0) {
		if (Outfq==1) {
			cerr<<"Fasta format input file can't output fastq format file"<<endl;
			return 1;
		}
		Run_fasta(P2In, Infq, Outgz);
	}
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

