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
#include <cstring> // for strlen
#include <unordered_map>
#include <cassert>
#include <map>
#include "zlib.h"
#include "kseq.h"
#include "kc-c4.c"
#include <cstdio>
#include "edlib.cpp"
#include "libdeflate.h"

#define OUTPUT_BUFFER_SIZE  1048576

using namespace std;

//KSEQ_INIT(gzFile, gzread)

int  TGSFilter_usage() {
	cout <<""
		"Usage: tgsfilter -1 TGS.raw.fq.gz -o TGS.clean.fq.gz\n"
		" Input/Output options:\n"
		"   -i	<str>   input of fasta/q file\n"
		"   -o	<str>   output of fasta/q file\n"
		" Basic filter options:\n"
		"   -l	<int>   min length of read to out [1000]\n"
		"   -L	<int>   max length of read to out\n"
		"   -q	<int>   min mean base quality [auto]\n"
		"   -5	<int>   trim bases from the 5' end of the read [0]\n"
		"   -3	<int>   trim bases from the 3' end of the read [0]\n"
		" Adapter filter options:\n"
		"   -a	<str>   adapter sequence file \n"
		"   -A           disable reads filter, only for adapter identify\n"
		"   -N	<int>   read number for adapter identify [200000]\n"
		"   -E	<int>   read end length for adapter identify [100]\n"
		"   -k	<int>   kmer size for adapter assembly [19]\n"
		"   -y	<int>   min assembly adapter length [30]\n"
		"   -m	<int>   min match length between read end and adapter [4]\n"
		"   -M	<int>   min match length between read middle and adapter [35]\n"
		"   -s  <float>  min similarity between read end and adapter [0.75]\n"
		"   -S  <float>  min similarity between read middle and adapter [0.9]\n"
		"   -D           split reads with middle adapter instead of discard\n"
		" Other options:\n"
		"   -t           number of threads [3]\n"
		"   -h           show help [v1.08]\n"
		"\n";
	return 1;
}

int n_thread=3;
uint64_t BLOCK_SIZE=512000; // bp

class Para_A24 {
	public:
		//
		string InFile;
		string OutFile;
		//
		int MinLength;
		uint64_t MaxLength;
		int AverQ;
		int HeadCrop;
		int TailCrop;
		//
		string AdapterFile;
		int ReadNumber;
		int EndLen;
		int Kmer;
		int AdapterLen;
		int EndMatchLen;
		int MidMatchLen;
		float EndSim;
		float MidSim;
		bool discard;
		//
		int Infq;
		int Outfq;
		bool OUTGZ;
		bool ONLYAD;
		int ReadLength;

		Para_A24() {
			InFile="";
			OutFile="";
			
			MinLength=1000;
			MaxLength=UINT64_MAX;
			AverQ=-1;
			HeadCrop=0;
			TailCrop=0;
			
			AdapterFile="";
			ReadNumber=200000;
			EndLen=100;
			Kmer=19;
			AdapterLen=30;
			EndMatchLen=4;
			MidMatchLen=35;
			EndSim=0.75;
			MidSim=0.9;
			discard=true;
			
			Infq=1;
			Outfq=1;
			OUTGZ=false;
			ONLYAD=false;
			ReadLength=0;
		}
};

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
			P2In->MinLength=atoi(argv[i]);
			if (P2In->MinLength<100){
				P2In->MinLength=100;
			}
		}
		else if (flag == "L") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MaxLength=atoi(argv[i]);
		}
		else if (flag == "q") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->AverQ=atoi(argv[i]);
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
		else if (flag == "k"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->Kmer=atoi(argv[i]);
		}
		else if (flag == "y"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->AdapterLen=atoi(argv[i]);
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
		else if (flag == "s"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndSim=atof(argv[i]);
			if (P2In->EndSim<0.7){
				P2In->EndSim=0.7;
				cout << "Warning: re set -s to : "<<P2In->EndSim<<endl;
			}
		}
		else if (flag == "S"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MidSim=atof(argv[i]);
			if (P2In->MidSim<0.85){
				P2In->MidSim=0.85;
				cout << "Warning: re set -S to : "<<P2In->MidSim<<endl;
			}
		}
		else if (flag == "D"){
			P2In->discard=false;
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
			BLOCK_SIZE=atoi(argv[i]);
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
	if (!(P2In->ONLYAD)){
		if ( (P2In->OutFile).empty()) {
			cerr<< "Error: -o lack argument for the must"<<endl;
			return 1;
		}
	}

	if (access((P2In->InFile).c_str(), 0) != 0) {
		cerr<<"Error: Can't find this file for -i "<<(P2In->InFile)<<endl;
		return 1;
	}
	else {
		return  0;
	}
}

class DeflateCompress{
	public:
		DeflateCompress( ) {
			m_compressor = libdeflate_alloc_compressor(6);
		}

		~DeflateCompress() {
			libdeflate_free_compressor(m_compressor);
		}

		bool compressData(const void* input, size_t inputSize, uint8_t  ** m_outputBuffer, size_t & compressedSize, size_t & OUT_BUFFER_SIZE,int &  Thread) 
		{
			if (OUT_BUFFER_SIZE<0.8*inputSize) { // compression >=80%
				OUT_BUFFER_SIZE=inputSize;
				delete [] m_outputBuffer[Thread];
				m_outputBuffer[Thread]=new uint8_t[OUT_BUFFER_SIZE];
			}

			compressedSize = libdeflate_gzip_compress(m_compressor, input, inputSize, m_outputBuffer[Thread], OUT_BUFFER_SIZE);
			if (compressedSize == 0) {
				return false;
			}
			else {
				return true;
			}
		}

	private:
		libdeflate_compressor* m_compressor;
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

char complement[256];
string rev_comp_seq(const string& dna) {
	
	string reverse_complement;
	for (int i = dna.size() - 1; i >= 0; --i) {
		reverse_complement += complement[dna[i]];
	}
	return reverse_complement;
}

string adapterSearch(Para_A24 *P2In, string &readsName, float &meanDep){
	/*
	PB-1 RSII/Sequel/Revio
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

	float sim=0.9;
	gzFile fp = gzopen(readsName.c_str(), "r");
  	kseq_t *ks = kseq_init(fp);

	while (kseq_read(ks) >= 0) {
		string ts=ks->seq.s;
        int tsLen=ks->seq.l;
		string bestName;
		string bestStrand;
		int minLen=0;
        for (const auto& pair : adapterLib) {
            string name = pair.first;
			string qsFor=pair.second.first;
            string qsRev=pair.second.second;
            int qsForLen=qsFor.length();
            int qsRevLen=qsRev.length();
            int minK=static_cast<int>((1-sim)*qsForLen)+1;
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
	kseq_destroy(ks);
  	gzclose(fp);

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
		if (meanDep<3){
			adapter="";
			strand="";
		}
    }

	if (strand=="-"){
		adapter=rev_comp_seq(adapter);
	}
	
	return adapter;
}

static int IntKmer2Value(const kc_c4x_t *h , uint64_t x) 
{
	int mask = (1<<(h->p)) - 1;
	kc_c4_t *g=h->h[x&mask];
	khint_t k;
	int absent;
	int count=0;
	uint64_t key = (x >> (h->p)) << KC_BITS;
	k = kc_c4_put(g, key, &absent);
	if (kh_exist(g, k)) {
		count = kh_key(g, k) & KC_MAX;
	}
	return count;
}

bool checkKmer(const std::string& kmer, int &count, int &minCount, int &lastCount, 
				std::unordered_set<std::string> &seen_kmers){
    if (count >= minCount) {
        int GC_count = std::count(kmer.begin(), kmer.end(), 'G') + std::count(kmer.begin(), kmer.end(), 'C');
        float ratio = static_cast<float>(GC_count) / kmer.length();
        if (ratio >= 0.2 && ratio <= 0.8) {
			if (seen_kmers.find(kmer) == seen_kmers.end()){
				if (lastCount==0){
            		return true;
				}else{
					int minCount = count < lastCount? count : lastCount;
					int maxCount = count > lastCount? count : lastCount;
					float ratio = static_cast<float>(maxCount) / minCount;
					if (ratio<4){
						return true;
					}

				}
			}
        }
    }
    return false;
}

bool checkCanAdapter(const std::string& adapter, std::unordered_set<std::string> &seen_adapters) {
    int aCount = std::count(adapter.begin(), adapter.end(), 'A');
    int gCount = std::count(adapter.begin(), adapter.end(), 'G');
    int cCount = std::count(adapter.begin(), adapter.end(), 'C');
    int tCount = std::count(adapter.begin(), adapter.end(), 'T');

	string adapter_rev=rev_comp_seq(adapter);
	bool flagAdapter_for=seen_adapters.find(adapter) == seen_adapters.end();
	bool flagAdapter_rev=seen_adapters.find(adapter_rev) == seen_adapters.end();

    if (aCount >= 2 && gCount >= 2 && cCount >= 2 && tCount >= 2 && flagAdapter_for && flagAdapter_rev) {
        return true;
    }
    
    return false;
}

string AssemblyAdapter(Para_A24 *P2In, const kc_c4x_t *h, string &FilePath, float &meanDep,
						std::vector<std::pair<std::string, int>> &canAdapters){
	gzFile fp = gzopen(FilePath.c_str(), "r");
  	kseq_t *ks = kseq_init(fp);

	int minCount=(P2In->ReadNumber)/100000;
	if (minCount<5){
		minCount=5;
	}

    int minLen=(P2In->AdapterLen);
    std::vector<std::pair<std::string, int>> canSeqs;
	std::unordered_set<std::string> seen_seqs;

	int i, l;
    int k=P2In->Kmer;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;

	while(kseq_read(ks) >= 0){
		string name=ks->name.s;
		string seq=ks->seq.s;
        int seqLen = ks->seq.l;
		std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        int adapterDep=0;
		int lastCount=0;
		string adapter;
		std::unordered_set<std::string> seen_kmers;
		
		for (i = l = 0, x[0] = x[1] = 0; i < seqLen; ++i) {
			int c = seq_nt4_table[(uint8_t)seq[i]];
			if (c < 4) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					string kmer=seq.substr(i-k+1,k);
					uint64_t y = x[0] < x[1]? x[0] : x[1];
					uint64_t hash_value = hash64(y, mask);
					int Count=IntKmer2Value(h,hash_value);
					bool flagKmer=checkKmer(kmer, Count, minCount, lastCount, seen_kmers);
					
					if (flagKmer){
						if (adapter.length()==0){
							adapter=kmer;
						}else{
							adapter+=kmer.substr(k-1);
						}
						adapterDep+=Count;
						seen_kmers.insert(kmer);
						lastCount=Count;
					}else{
						bool flagAdapter=checkCanAdapter(adapter, seen_seqs);
						if (adapter.length()>=minLen && flagAdapter){
							canSeqs.push_back({adapter, adapterDep});
							seen_seqs.insert(adapter);
						}
						adapter="";
						adapterDep=0;
						lastCount=0;
						seen_kmers.clear();
					}
				}
			} 
			else {
				l = 0, x[0] = x[1] = 0;
			}
		}

		bool flagAdapter=checkCanAdapter(adapter, seen_seqs);
        if (adapter.length()>=minLen && flagAdapter){
            canSeqs.push_back({adapter, adapterDep});
			seen_seqs.insert(adapter);
        }
	}
	kseq_destroy(ks);
  	gzclose(fp);

	//
	std::sort(canSeqs.begin(), canSeqs.end(),
              [](const std::pair<std::string, int> &a, const std::pair<std::string, int> &b) {
                  return a.second > b.second;
              });
	//
    //std::vector<std::pair<std::string, int>> canAdapters;
	std::unordered_set<std::string> seen_adapters;

	for (int j = 0; j < canSeqs.size(); j++){
		string seq=canSeqs[j].first;
		int depth=canSeqs[j].second;
		int seqLen=seq.length();
		int minKmerCount=int(depth/((seqLen-k+1)*10));
		
		int adapterDep=0;
		string adapter;

		for (i = l = 0, x[0] = x[1] = 0; i < seqLen; ++i) {
			int c = seq_nt4_table[(uint8_t)seq[i]];
			if (c < 4) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					string kmer=seq.substr(i-k+1,k);
					uint64_t y = x[0] < x[1]? x[0] : x[1];
					uint64_t hash_value = hash64(y, mask);
					int Count=IntKmer2Value(h,hash_value);
					
					if (Count>=minKmerCount){
						if (adapter.length()==0){
							adapter=kmer;
						}else{
							adapter+=kmer.substr(k-1);
						}
						adapterDep+=Count;
					}else{
						bool flagAdapter=checkCanAdapter(adapter, seen_adapters);
						if (adapter.length()>=minLen && flagAdapter){
							canAdapters.push_back({adapter, adapterDep});
							seen_adapters.insert(adapter);
						}
						adapter="";
						adapterDep=0;
					}
				}
			} 
			else {
				l = 0, x[0] = x[1] = 0;
			}
		}

		bool flagAdapter=checkCanAdapter(adapter, seen_adapters);
        if (adapter.length()>=minLen && flagAdapter){
            canAdapters.push_back({adapter, adapterDep});
			seen_adapters.insert(adapter);
        }
	}
    //
	std::sort(canAdapters.begin(), canAdapters.end(),
              [](const std::pair<std::string, int> &a, const std::pair<std::string, int> &b) {
                  return a.second > b.second;
              });
    
	string best_adapter = canAdapters[0].first;
	meanDep = static_cast<float>(canAdapters[0].second) / (best_adapter.length() - k + 1);

	return best_adapter;
}

string checkAdapter(string &bestAdapter, int &totalDep, std::vector<std::pair<std::string, int>> &canAdapters){
	int checkSize=1000;
	if (canAdapters.size()<checkSize){
		checkSize=canAdapters.size();
	}

	string tsFor=bestAdapter;
	string tsRev=rev_comp_seq(tsFor);
	int tsLen=tsFor.length();

	int minLen=0;
	int minDep=0;
	string adapter;

	for (int j = 0; j < checkSize; j++){
		int depth=canAdapters[j].second;
		string qs=canAdapters[j].first;
		int qsLen=qs.length();
		EdlibAlignResult result_for = edlibAlign(qs.c_str(), qsLen, tsFor.c_str(), tsLen, 
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result_for.status == EDLIB_STATUS_OK){
			int dist=result_for.editDistance;
			int length=result_for.alignmentLength;
			int numAln=result_for.numLocations;
			int ts=result_for.startLocations[0];
			int te=result_for.endLocations[0]+1;
			int mlen=length-dist;
			float sim=static_cast<float>(mlen) / (te-ts);
			if(sim>=0.9){
				if (mlen > minLen){
					minLen=mlen;
					adapter=qs;
					minDep=depth;
				}
			}
		}
		edlibFreeAlignResult(result_for);
		//
		EdlibAlignResult result_rev = edlibAlign(qs.c_str(), qsLen, tsRev.c_str(), tsLen, 
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result_rev.status == EDLIB_STATUS_OK){
			int dist=result_rev.editDistance;
			int length=result_rev.alignmentLength;
			int numAln=result_rev.numLocations;
			int ts=result_rev.startLocations[0];
			int te=result_rev.endLocations[0]+1;
			int mlen=length-dist;
			float sim=static_cast<float>(mlen) / (te-ts);
			if(sim>=0.9){
				if (mlen > minLen){
					minLen=mlen;
					adapter=qs;
					minDep=depth;
				}
			}
		}
		edlibFreeAlignResult(result_rev);
	}

	string outAdapter;
	if (minLen > (adapter.length()/2)){
		totalDep=minDep;
		outAdapter=adapter;
	}
	return outAdapter;
}

std::map<std::string, std::pair<std::string, std::string>> adapters;
bool asseFlag=false;

int Get_filter_parameter(Para_A24 *P2In){
	string FilePath;
	string InPath=(P2In->InFile);
	string OutPath=(P2In->OutFile);
	if (OutPath.empty()){
		FilePath=InPath;
	}else{
		FilePath=OutPath;
	}

	gzFile fp;
  	kseq_t *seq;

	fp = gzopen(InPath.c_str(), "r");
  	seq = kseq_init(fp);

	int seqNum=0;
	int maxSeq=(P2In->ReadNumber);
	
	int qtNum=10000;
	int minQ=50000;
	int maxQ=0;
	uint64_t sumQ = 0;
	uint64_t sumLen = 0;

	std::map<std::string, std::string> headSeq;
	std::map<std::string, std::string> tailSeq;

	while((kseq_read(seq)) >= 0){
		if (seqNum>=maxSeq){
			break;
		}

		uint64_t length=seq->seq.l;
		if (length<5000){
			continue;
		}
		seqNum++;
		//
		if (seqNum < qtNum){
			string qual=seq->qual.s;
			uint64_t qual_len=seq->qual.l;
			sumLen+=qual_len;
			for (int i=0; i<qual_len; i++) {
				sumQ+=qual[i];
				if(minQ>qual[i]) {
					minQ=qual[i];
				}
				if(maxQ<qual[i]) {
					maxQ=qual[i];
				}
			}
		}

		string name=seq->name.s;
		string sequence=seq->seq.s;

		string head_seq=sequence.substr(0, P2In->EndLen);
		string tail_seq=sequence.substr(length-(P2In->EndLen));
	
		headSeq[name] = head_seq;
		tailSeq[name] = tail_seq;
	}

	kseq_destroy(seq);
  	gzclose(fp);

	//mean phred quality
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
	if(P2In->AverQ<0){
		if (meanQ>=25){
			P2In->AverQ=20;
		}else if(meanQ>0) {
			P2In->AverQ=10;
		}else{
			P2In->AverQ=0;
		}
	}

	//search adapter
	if (!(P2In->AdapterFile).empty() && access((P2In->AdapterFile).c_str(), 0) == 0){
		string adapter_name=P2In->AdapterFile;
		gzFile fp = gzopen(adapter_name.c_str(), "r");
    	kseq_t *ks = kseq_init(fp);
		while (kseq_read(ks) >= 0) {
			string name=ks->name.s;
			string seq_for=ks->seq.s;
			string seq_rev=rev_comp_seq(seq_for);
			adapters[name] = {seq_for,seq_rev};
		}
		kseq_destroy(ks);
		gzclose(fp);
		return 0;
	}

	string prefix=GetFilePreifx(FilePath);
	string out_front=prefix + ".front.fa";
	ofstream OUTHF(out_front);
	for (const auto &info : headSeq) {
		OUTHF<<">"<<info.first<<"\n"<<info.second<<endl;
	}
	OUTHF.close();

	string out_tail=prefix + ".tail.fa";
	ofstream OUTHT(out_tail);
	for (const auto &info : tailSeq) {
		OUTHT<<">"<<info.first<<"\n"<<info.second<<endl;
	}
	OUTHT.close();

	headSeq.clear();
	tailSeq.clear();
	std::map<std::string, std::string>().swap(headSeq);
	std::map<std::string, std::string>().swap(tailSeq);

	string adapter_5p;
	string adapter_3p;
	float depth_5p=0;
	float depth_3p=0;
	int len_5p=0;
	int len_3p=0;
	std::vector<std::pair<std::string, int>> can5pAdapters;
	std::vector<std::pair<std::string, int>> can3pAdapters;

	cout << "INFO: searching 5' adapter..."<<endl;
	adapter_5p=adapterSearch(P2In, out_front, depth_5p);
	len_5p=adapter_5p.length();
	
	cout << "INFO: searching 3' adapter..."<<endl;
	adapter_3p=adapterSearch(P2In, out_tail, depth_3p);
	len_3p=adapter_3p.length();

	if(len_5p==0 && len_3p==0){
		asseFlag=true;
		//cout << "Warnings: the adapter could not be found within the general adapter library."<<endl;
		int p=15;
		uint64_t block_size = 10000000;

		cout << "INFO: assembly 5' adapter..."<<endl;
		kc_c4x_t *fh;
		fh = count_file(out_front.c_str(), P2In->Kmer, p, block_size, n_thread);
		adapter_5p=AssemblyAdapter(P2In, fh, out_front, depth_5p, can5pAdapters);
		len_5p=adapter_5p.length();
		for (int i = 0; i < 1<<p; ++i) kc_c4_destroy(fh->h[i]);
		free(fh->h); free(fh);

		cout << "INFO: assembly 3' adapter..."<<endl;
		kc_c4x_t *th;
		th = count_file(out_tail.c_str(), P2In->Kmer, p, block_size, n_thread);
		adapter_3p=AssemblyAdapter(P2In, th, out_tail, depth_3p, can3pAdapters);
		len_3p=adapter_3p.length();
		for (int i = 0; i < 1<<p; ++i) kc_c4_destroy(th->h[i]);
		free(th->h); free(th);
	}

	//check depth 
	if (depth_5p > 5*depth_3p){
		adapter_3p="";
		depth_3p=0;
	}else if(depth_3p > 5*depth_5p){
		adapter_5p="";
		depth_5p=0;
	}else{
		string outAdapter;
		int outLen;
		int totalDep;
		if (depth_5p > depth_3p){
			outAdapter=checkAdapter(adapter_5p, totalDep, can3pAdapters);
			outLen=outAdapter.length();
			if(outLen > 0){
				adapter_3p=outAdapter;
				depth_3p=static_cast<float>(totalDep)/(outLen+1-(P2In->Kmer));
			}
		}else if(depth_5p < depth_3p){
			outAdapter=checkAdapter(adapter_3p, totalDep, can5pAdapters);
			outLen=outAdapter.length();
			if(outLen > 0){
				adapter_5p=outAdapter;
				depth_5p=static_cast<float>(totalDep)/(outLen+1-(P2In->Kmer));
			}
		}
	}

	cout <<"INFO: 5' adapter: " << adapter_5p << endl;
	cout <<"INFO: 3' adapter: " << adapter_3p << endl;

	cout <<"INFO: mean depth of 5' adapter: "<<depth_5p<<endl;
	cout <<"INFO: mean depth of 3' adapter: "<<depth_3p<<endl;

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
	
	remove(out_front.c_str());
	remove(out_tail.c_str());
	//
	
	if (!(P2In->ONLYAD)){
		cout << "INFO: trim 5' end length: "<<P2In->HeadCrop<<endl;
		cout << "INFO: trim 3' end length: "<<P2In->TailCrop<<endl;
		if ((P2In->Infq)==1){
			cout<<"INFO: min mean base quality was set to: "<<(P2In->AverQ)<<endl;
			P2In->AverQ=(P2In->AverQ)+qType;
		}
		cout << "INFO: min output reads length: "<<P2In->MinLength<<endl;
	}
	

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

int GetEditDistance(Para_A24 *P2In, string &query, string &target, std::vector<std::vector<int>> &adapterRegions){
	
	int midNum=0;
	float endSim=P2In->EndSim;
	float midSim=P2In->MidSim;
	int qLen=query.length();
	int tLen=target.length();

	int end5Len=(P2In->EndLen) - (P2In->HeadCrop);
	int end3Len=(P2In->EndLen) - (P2In->TailCrop);

	int maxK=int((1-midSim)*qLen)+1;

	if (asseFlag){
		maxK=qLen-(P2In->MidMatchLen)+1;
	}

	//check middle adapter
	EdlibAlignResult result = edlibAlign(query.c_str(), qLen, target.c_str(), tLen, 
                    edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	if (result.status == EDLIB_STATUS_OK){
		int dist=result.editDistance;
		int numAln=result.numLocations;
		int length=result.alignmentLength;
		int mlen=length - dist;
		if (mlen >= (P2In->MidMatchLen)){
			for (int i=0; i<numAln; i++){
				int ts=result.startLocations[i];
				int te=result.endLocations[i]+1;
				float sim = static_cast<float>(mlen) / length;
				if (sim >= midSim && ts>=end5Len && (tLen-te)>=end3Len){
					//cout <<tLen<<" "<<length<<" mid "<<dist<<" "<<mlen<<" "<<ts<<" "<<te<<endl;
					midNum++;
					adapterRegions.push_back({ts, te});
				}
			}
		}
	}
	edlibFreeAlignResult(result);

	//check 5' and 3' end
	int checkLen= (P2In->EndLen) + qLen;
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
						//cout <<tLen<<" "<<length<<" 5p "<<dist<<" "<<mlen<<" "<<ts<<" "<<te<<endl;
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
			for (int i=0; i<numAln; i++){
				int ts=result_ts3.startLocations[i] + tLen - ts3Len;
				int te=result_ts3.endLocations[i]+1 + tLen - ts3Len;
				if (mlen >= (P2In->EndMatchLen)){
					float sim1 = static_cast<float>(mlen) / (tLen-ts);
					float sim2 = static_cast<float>(mlen) / length;
					if (sim1 >= endSim || sim2 >= endSim){
						//cout <<tLen<<" "<<length<<" 3p "<<dist<<" "<<mlen<<" "<<ts<<" "<<te<<endl;
						adapterRegions.push_back({ts, tLen});
					}
				}
			}
		}
		edlibFreeAlignResult(result_ts3);
	}
	return midNum;
}

void adapterMap(Para_A24 * P2In, string &raw_seq, int &raw_len,
				std::vector<std::vector<int>> &keepRegions,
				std::array<uint64_t, 16>& DropInfo){

	std::vector<std::vector<int>> adapterRegions;
	int mid_for;
	int mid_rev;
	//
    for (const auto& pair : adapters) {
        string qsFor=pair.second.first;
        string qsRev=pair.second.second;
		mid_for = GetEditDistance(P2In, qsFor, raw_seq, adapterRegions);
		mid_rev = GetEditDistance(P2In, qsRev, raw_seq, adapterRegions);
    }

	int midNum=0;
	if (mid_for > 0 || mid_rev > 0){
		midNum=1;
	}
	DropInfo[14]+=midNum;

	if (midNum>0 && P2In->discard){
		DropInfo[15]+=raw_len;
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
					if (keep_len >= (P2In->MinLength) && keep_len <= (P2In->MaxLength)){
						keepRegions.push_back({currentStart, keep_len});
					}else{
						DropInfo[12]++;
						DropInfo[13]+=keep_len;
					}
				}
				currentStart = region[1];
			}

			if (currentStart < raw_len) {
				keep_len=static_cast<int>(raw_len) - currentStart;
				if (keep_len >= (P2In->MinLength) && keep_len <= (P2In->MaxLength)){
					keepRegions.push_back({currentStart, keep_len});
				}else{
					DropInfo[12]++;
					DropInfo[13]+=keep_len;
				}
			}
		}else{
			if (raw_len >= (P2In->MinLength) && raw_len <= (P2In->MaxLength)){
				keepRegions.push_back({currentStart, raw_len});
			}else{
				DropInfo[12]++;
				DropInfo[13]+=raw_len;
			}
		}
	}
}

void Filter_fastq_reads_adapter(Para_A24 * P2In, string &OUT_DATA,
						vector<string>& ID, vector <string> &SEQ, vector <string> &QUAL,
						std::array<uint64_t, 16>& DropInfo){

	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);
	
	for (int i = 0; i < SEQ.size(); i++) {
		int seq_len = SEQ[i].length();
		int qual_len = QUAL[i].length();
		DropInfo[0]++;
		DropInfo[1]+=seq_len;

		if (seq_len < (P2In->MinLength + totalCrop) || (seq_len != qual_len)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
			continue;
		}
		
		string raw_seq = SEQ[i].substr(headCrop, seq_len - totalCrop);
		string raw_qual = QUAL[i].substr(headCrop, qual_len - totalCrop);
		int raw_len = raw_seq.length();

		if (totalCrop>0){
			DropInfo[6]++;
			DropInfo[7]+=totalCrop;
		}

		//
		double avgQuality = CalcAvgQuality(raw_qual);
		if (avgQuality < (P2In->AverQ)){
			DropInfo[10]++;
			DropInfo[11]+=raw_len;
			continue;
		}
		
		std::vector<std::vector<int>> keepRegions;
		adapterMap(P2In, raw_seq, raw_len, keepRegions,DropInfo);
		
		string raw_id=ID[i];
		string seq_id=raw_id;
		int pass_num=0;
		if ((P2In->Outfq)==1){
			if (keepRegions.size()>0){
				for (const auto& region : keepRegions) {
					pass_num++;
					int start = region[0];
					int seqLen = region[1];
					string seq=raw_seq.substr(start,seqLen);
					string qual=raw_qual.substr(start,seqLen);

					DropInfo[2]++;
					DropInfo[3]+=seqLen;

					if (pass_num>=2){
						seq_id=raw_id+":"+std::to_string(pass_num);
					}
					OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + seq + 
							 "\n+\n" + qual + "\n";
				}
			}
		}else{
			if (keepRegions.size()>0){
				for (const auto& region : keepRegions) {
					pass_num++;
					int start = region[0];
					int seqLen = region[1];
					string seq=raw_seq.substr(start,seqLen);
					string qual=raw_qual.substr(start,seqLen);

					DropInfo[2]++;
					DropInfo[3]+=seqLen;

					if (pass_num>=2){
						seq_id=raw_id+":"+std::to_string(pass_num);
					}
					OUT_DATA=OUT_DATA+">"+seq_id+"\n" + seq + "\n";
				}
			}
		}
	}
}

void Filter_fasta_reads_adapter(Para_A24 * P2In,string &OUT_DATA, 
						vector<string>& ID, vector <string> & SEQ,
						std::array<uint64_t, 16>& DropInfo) {
	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);
	
	for (int i = 0; i < SEQ.size(); i++) {
		int seq_len = SEQ[i].length();
		DropInfo[0]++;
		DropInfo[1]+=seq_len;

        if (seq_len < (P2In->MinLength + totalCrop)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
            continue;
        }

        string raw_seq = SEQ[i].substr(headCrop, seq_len - totalCrop);
        int raw_len = raw_seq.length();

		if (totalCrop>0){
			DropInfo[6]++;
			DropInfo[7]+=totalCrop;
		}

        std::vector<std::vector<int>> keepRegions;
		adapterMap(P2In, raw_seq, raw_len, keepRegions,DropInfo);
        
        string seq_id=ID[i];
        int pass_num=0;
		if (keepRegions.size()>0){
			for (const auto& region : keepRegions) {
				pass_num++;
				int start = region[0];
				int seqLen = region[1];
				string seq=raw_seq.substr(start,seqLen);
		
				DropInfo[2]++;
				DropInfo[3]+=seqLen;

				if (pass_num>=2){
					seq_id=ID[i]+":"+std::to_string(pass_num);
				}
				
				OUT_DATA=OUT_DATA+">"+seq_id+"\n" + seq + "\n";
			}
		}
	}
}

void Filter_reads_adapter(Para_A24 * P2In, string &OUT_DATA, 
				vector<string> &ID, vector <string> &SEQ, vector <string> &QUAL,
				std::array<uint64_t, 16>& DropInfo){
	OUT_DATA="";
	if(QUAL[0].length()<=2){
		Filter_fasta_reads_adapter(P2In, OUT_DATA, ID, SEQ, DropInfo);
	} else {
		Filter_fastq_reads_adapter(P2In, OUT_DATA, ID, SEQ, QUAL, DropInfo);
	}
}

void Filter_reads_adapter_gz(Para_A24 *P2In, vector<string> &ID, vector <string> &SEQ, 
					vector <string> &QUAL,
					std::array<uint64_t, 16>& DropInfo,
					uint8_t ** ComData, size_t & ComSize,
					size_t &ComBuff, int &tid){
	string OUT_DATA="";
	OUT_DATA.reserve(OUTPUT_BUFFER_SIZE);
	if(QUAL[0].length()<=2){
		Filter_fasta_reads_adapter(P2In, OUT_DATA, ID, SEQ, DropInfo);

	} else {
		Filter_fastq_reads_adapter(P2In, OUT_DATA, ID, SEQ, QUAL, DropInfo);
	}
	
	if (!(OUT_DATA.empty())) {
		DeflateCompress  GZData;
		size_t inputSize  = OUT_DATA.length();
		GZData.compressData(OUT_DATA.c_str(), inputSize, ComData, 
							ComSize,ComBuff,tid);
	}else{
		ComSize=0;
	}
}

int Run_seq_filter_adapter(Para_A24 * P2In) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	std::vector<std::vector<std::string>> ID(n_thread);
	std::vector<std::vector<std::string>> SEQ(n_thread);
	std::vector<std::vector<std::string>> QUAL(n_thread);

	std::vector<std::array<uint64_t, 16>> DropInfo(n_thread, std::array<uint64_t, 16>{});
	//raw_reads raw_bases
	//clean_reads clean_bases
	//shortDrop_reads shortDrop_bases
	//endDrop_reads  endDrop_bases
	//adapterDrop_reads adapterDrop_bases
	//LQDrop_reads LQDrop_bases
	//OutDrop_reads OutDrop_bases
	//Middle_adapter drop_bases

	std::vector<std::thread> threads;

	gzFile fp = gzopen((P2In->InFile).c_str(), "r");
    kseq_t *seq = kseq_init(fp);

	string outname=P2In->OutFile;
	int total_length = 0;
	//

	if (P2In->OUTGZ){
		ofstream OUTHGZ(outname.c_str(), std::ios::out | std::ios::binary);
		uint8_t ** ComData = new uint8_t*[n_thread];
		size_t * ComSize =new size_t [n_thread];
		size_t * ComBuff =new size_t [n_thread];
		
		int * ArryThread=new int [n_thread];
		for (int i = 0; i < n_thread; i++) {
			ComData[i] = new uint8_t[OUTPUT_BUFFER_SIZE];
			ComSize[i]=OUTPUT_BUFFER_SIZE;
			ComBuff[i]=OUTPUT_BUFFER_SIZE;
			ArryThread[i]=i;
		}
		
		int tid=0;

		while (kseq_read(seq) >= 0) {
			string name=seq->name.s;
			string sequence=seq->seq.s;
			string qual=seq->qual.s;
			int seq_length=seq->seq.l;
			
			total_length += seq_length;

			ID[tid].push_back(name);
			SEQ[tid].push_back(sequence);
			QUAL[tid].push_back(qual);
			
			if (total_length >= BLOCK_SIZE) {
				threads.push_back(thread(Filter_reads_adapter_gz,P2In, ref(ID[tid]), 
									ref(SEQ[tid]),ref(QUAL[tid]),
									ref(DropInfo[tid]), ComData, ref(ComSize[tid]),
									ref(ComBuff[tid]), ref(ArryThread[tid])));
				tid++;
				total_length = 0;
				
				if (threads.size() == n_thread) {
					for (auto &t : threads) {
						t.join();
					}
					threads.clear();
					
					for (int i = 0; i < n_thread; i++) {
						if (ComSize[i]>0) {
							OUTHGZ.write((const char*) ComData[i], ComSize[i] );
						}
						ID[i].clear();
						SEQ[i].clear();
						QUAL[i].clear();
					}
					tid=0;
				}
			}	
		}

		if (!SEQ[tid].empty()) {
			threads.push_back(thread(Filter_reads_adapter_gz,P2In, ref(ID[tid]), 
							ref(SEQ[tid]),ref(QUAL[tid]),
							ref(DropInfo[tid]), ComData, ref(ComSize[tid]),
							ref(ComBuff[tid]), ref(ArryThread[tid])));
			total_length = 0;
		}
		
		for (auto &t : threads) {
			t.join();
		}

		threads.clear();

		for (int i = 0; i <= tid; i++) {
			if (ComSize[i]>0) {
				OUTHGZ.write((const char*) ComData[i], ComSize[i] );
			}
			ID[i].clear();
			SEQ[i].clear();
			QUAL[i].clear();
		}
		
		for (int i = 0; i < n_thread; i++) {
			delete[]  ComData[i] ;
		}
			
		delete [] ComSize;
		delete [] ComBuff;
		delete [] ComData;
		delete [] ArryThread;

	} else {
		ofstream OUTH;
		OUTH.open(outname.c_str());
		vector <string> OUT_DATA;
		OUT_DATA.resize(n_thread);

		int tid=0;

		while (kseq_read(seq) >= 0) {
			string name=seq->name.s;
			string sequence=seq->seq.s;
			string qual=seq->qual.s;
			int seq_length=seq->seq.l;

			total_length += seq_length;

			ID[tid].push_back(name);
			SEQ[tid].push_back(sequence);
			QUAL[tid].push_back(qual);

			if (total_length >= BLOCK_SIZE) {
				
				threads.push_back(thread(Filter_reads_adapter,P2In, ref(OUT_DATA[tid]), 
						ref(ID[tid]),ref(SEQ[tid]),ref(QUAL[tid]), ref(DropInfo[tid])));
				tid++;
				total_length = 0;

				if (threads.size() == n_thread) {
					for (auto &t : threads) {
						t.join();
					}
					threads.clear();
					for (int i = 0; i < n_thread; i++) {
						if (!(OUT_DATA[i].empty())) {
							OUTH << OUT_DATA[i];
						}
						ID[i].clear();
						SEQ[i].clear();
						QUAL[i].clear();
					}
					tid=0;
				}
			}
		}

		if (!SEQ[tid].empty()) {
			threads.push_back(thread(Filter_reads_adapter,P2In, ref(OUT_DATA[tid]), 
					ref(ID[tid]),ref(SEQ[tid]),ref(QUAL[tid]), ref(DropInfo[tid])));
			total_length = 0;
		}

		for (auto &t : threads) {
			t.join();
		}

		threads.clear();

		for (int i = 0; i <= tid; i++) {
			if (!(OUT_DATA[i].empty())) {
				OUTH << OUT_DATA[i];
			}
			ID[i].clear();
			SEQ[i].clear();
			QUAL[i].clear();
		}
		
		
		OUTH.close();
	}
	kseq_destroy(seq);
	gzclose(fp);

	uint64_t raw_reads=0;
	uint64_t raw_bases=0;
	uint64_t clean_reads=0;
	uint64_t clean_bases=0;
	uint64_t shortDrop_reads=0;
	uint64_t shortDrop_bases=0;
	uint64_t endDrop_reads=0;
	uint64_t endDrop_bases=0;
	uint64_t adapterDrop_reads=0;
	uint64_t adapterDrop_bases=0;
	uint64_t LQDrop_reads=0;
	uint64_t LQDrop_bases=0;
	uint64_t outDrop_reads=0;
	uint64_t outDrop_bases=0;
	uint64_t middle_reads=0;
	uint64_t middle_bases=0;

	for (int i = 0; i < n_thread; i++) {
		raw_reads+=DropInfo[i][0];
		raw_bases+=DropInfo[i][1];
		clean_reads+=DropInfo[i][2];
		clean_bases+=DropInfo[i][3];
		shortDrop_reads+=DropInfo[i][4];
		shortDrop_bases+=DropInfo[i][5];
		endDrop_reads+=DropInfo[i][6];
		endDrop_bases+=DropInfo[i][7];
		adapterDrop_reads+=DropInfo[i][8];
		adapterDrop_bases+=DropInfo[i][9];
		LQDrop_reads+=DropInfo[i][10];
		LQDrop_bases+=DropInfo[i][11];
		outDrop_reads+=DropInfo[i][12];
		outDrop_bases+=DropInfo[i][13];
		middle_reads+=DropInfo[i][14];
		middle_bases+=DropInfo[i][15];
	}

	cout << "INFO: "<< raw_reads<<" reads with a total of "<<raw_bases<<" bases were input."<<endl;
	cout << "INFO: "<< shortDrop_reads <<" reads were discarded with "<<shortDrop_bases<<" bases before filtering."<<endl;
	cout << "INFO: "<< endDrop_reads <<" reads were trimmed by "<<endDrop_bases<<" bases at the end."<<endl;
	cout << "INFO: "<< LQDrop_reads <<" reads were discarded with "<<LQDrop_bases<<" bases due to low quality."<<endl;
	cout << "INFO: "<< middle_reads <<" read was discarded with " <<middle_bases<<" bases due to a middle adapter."<<endl;
	cout << "INFO: "<< adapterDrop_reads <<" reads were trimmed by " <<adapterDrop_bases<<" bases with an adapter."<<endl;
	cout << "INFO: "<< outDrop_reads<<" reads were discarded with "<<outDrop_bases<<" bases before output."<<endl;
	cout << "INFO: "<< clean_reads <<" reads with a total of "<<clean_bases<<" bases were output"<<endl;

	return 0;
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

	string OutPath;
	if (!(P2In->ONLYAD)){
		OutPath=(P2In->OutFile);	
		P2In->Outfq = GetFileType(OutPath);
		string ext = GetFileExtension(OutPath);
		if (ext=="gz"){
			P2In->OUTGZ=true;
		}
	}

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
		cerr<<"Error: Fasta format input file can't output fastq format file"<<endl;
		return 1;
	}

	//
	Get_filter_parameter(P2In);

	if (P2In->ONLYAD){
		return 0;
	}
	
	Run_seq_filter_adapter(P2In);

	delete P2In ;
	return 0;
}

