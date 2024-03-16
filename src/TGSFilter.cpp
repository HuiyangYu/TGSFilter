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
#include "minimap.h"
#include "libdeflate.cpp"

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
		"   -w	<int>   windows size to cut off low quality region [50]\n"
		"   -n	<int>   read number for base content check [200000]\n"
		"   -e	<int>   read end length for base content check [100]\n"
		"   -5	<int>   drop bases from the front (5') of the read [auto]\n"
		"   -3	<int>   drop bases from the tail (3') of the read [auto]\n"
		//"   -x	<int>   min length to detect polyX in the read tail [10]\n"
		//"   -d          disable polyA trimming \n"
		" Adapter filter options:\n"
		"   -a	<str>   adapter sequence file \n"
		"   -A           disable reads filter, only for adapter identify\n"
		"   -N	<int>   read number for adapter identify [200000]\n"
		"   -E	<int>   read end length for adapter identify [100]\n"
		"   -k	<int>   kmer size for adapter assembly [15]\n"
		"   -y	<int>   min assembly adapter length [20]\n"
		"   -m	<int>   min match length between read end and adapter [5]\n"
		"   -M	<int>   min match length between read middle and adapter [auto]\n"
		"   -s  <float>  min similarity between read and adapter [0.8]\n"
		" Other options:\n"
		"   -t           number of threads [3]\n"
		"   -h           show help [v1.07]\n"
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
		int Window;
		int BCNum;
		int BCLen;
		int HeadCrop;
		int TailCrop;
		int PolyX;
		bool PolyA;
		
		//
		string AdapterFile;
		int ReadNumber;
		int EndLen;
		int Kmer;
		int AdapterLen;
		int MatchLen;
		int MidLen;
		float Similarity;
		//
		int Infq;
		int Outfq;
		bool OUTGZ;
		bool ONLYAD;
		int ReadLength;

		Para_A24() {
			//
			InFile="";
			OutFile="";
			//
			MinLength=1000;
			MaxLength=UINT64_MAX;
			AverQ=-1;
			Window=50;
			BCNum=200000;
			BCLen=100;
			HeadCrop=-1;
			TailCrop=-1;
			PolyX=10;
			PolyA=false;
			
			//
			AdapterFile="";
			ReadNumber=200000;
			EndLen=100;
			Kmer=15;
			AdapterLen=20;
			MatchLen=5;
			MidLen=-1;
			Similarity=0.8;
			//
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
			if (P2In->MinLength<1){
				P2In->MinLength=1;
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
		else if (flag == "w"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->Window=atoi(argv[i]);
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
		else if (flag == "x"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->PolyX=atoi(argv[i]);
		}
		else if (flag == "d"){
			P2In->PolyX=true;
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
			P2In->MatchLen=atoi(argv[i]);
		}
		else if (flag == "M"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MidLen=atoi(argv[i]);
		}
		else if (flag == "s"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->Similarity=atof(argv[i]);
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

void get_seqDict(string &libname, std::map<std::string, std::string> &seqDict){
	gzFile fp = gzopen(libname.c_str(), "r");
    kseq_t *ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		seqDict[ks->name.s] = ks->seq.s;
	}
	kseq_destroy(ks);
	gzclose(fp);
}

string adapterSearch(Para_A24 *P2In, string &libname, string &readsname, int &mean_depth){
	
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;

	mm_verbose = 2;
	mm_set_opt(0, &iopt, &mopt);
	iopt.k=15;
	iopt.w=1;
	mopt.flag |= MM_F_CIGAR;
	mopt.min_cnt = 1;
	mopt.min_chain_score =15;
	mopt.min_dp_max=15;

	gzFile f = gzopen(readsname.c_str(), "r");
    assert(f);
    kseq_t *ks = kseq_init(f);

	mm_idx_reader_t *r = mm_idx_reader_open(libname.c_str(), &iopt, 0);
    mm_idx_t *mi;

	int index_thread=1;
	std::map<std::pair<std::string, std::string>, int> maps;

	while ((mi = mm_idx_reader_read(r, index_thread)) != 0) {
        mm_mapopt_update(&mopt, mi);
        mm_tbuf_t *tbuf = mm_tbuf_init();
        gzrewind(f);
        kseq_rewind(ks);
        while (kseq_read(ks) >= 0) {
			mm_reg1_t *reg;
			int j, i, n_reg;

			reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0);
			for (j = 0; j < n_reg; ++j) {
				mm_reg1_t *r = &reg[j];
				int ql=ks->seq.l;
				int qs=r->qs;
				int qe=r->qe;
				int tl=mi->seq[r->rid].len;
				int ts=r->rs;
				int te=r->re;
				int mlen=r->mlen;
				int oh=tl * (1-(P2In->Similarity));
				string name=mi->seq[r->rid].name;
				string strand = r->rev ? "-" : "+";
				float similarity;
				if (qs<=oh){
					similarity=static_cast<float>(mlen) / qe;
				}else if((ql-qe)<=oh){
					similarity=static_cast<float>(mlen) / (ql-qs);
				}else{
					int blen=r->blen;
					if (blen<tl){
						blen=tl;
					}
					similarity=static_cast<float>(mlen) / blen;
				}

				if (similarity>=(P2In->Similarity)){
					maps[std::make_pair(name, strand)]+=mlen;
				}
			}
			free(reg);
		}
		mm_tbuf_destroy(tbuf);
        mm_idx_destroy(mi);
	}
	mm_idx_reader_close(r); // close the index reader
    kseq_destroy(ks);
	gzclose(f);
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
	std::map<std::string, std::string> seqDict;
	get_seqDict(libname,seqDict);

	if (!maps_vec.empty()) {
        const auto& maxItem = maps_vec.front();
		name=maxItem.first.first;
		adapter=seqDict[name];
		strand=maxItem.first.second;
		mean_depth=maxItem.second/(adapter.length());
		if (mean_depth<5){
			adapter="";
			strand="";
		}
    }

	if (strand=="-"){
		adapter=rev_comp_seq(adapter);
	}
	if (adapter.length()>0){
		cout << "INFO: found adapter"<<endl;
		cout <<">"<<name << "\n"<<adapter<<endl;
	}else{
		cout << "INFO: not found adapter"<<endl;
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

bool checkKmer(const std::string& kmer, int &count, int &minCount){
    if (count >= minCount) {
        int GC_count = std::count(kmer.begin(), kmer.end(), 'G') + std::count(kmer.begin(), kmer.end(), 'C');
        float ratio = static_cast<float>(GC_count) / kmer.length();
        if (ratio >= 0.2 && ratio <= 0.8) {
            return true;
        }
    }
    return false;
}

string AssemblyAdapter(Para_A24 *P2In, const kc_c4x_t *h, string &FilePath, int &mean_depth){

	gzFile fp = gzopen(FilePath.c_str(), "r");
  	kseq_t *seq = kseq_init(fp);

	int i, l;
    int k=P2In->Kmer;
	int seq_len;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	
	int minCount=(P2In->ReadNumber)/100000;
	if (minCount<5){
		minCount=5;
	}

	std::vector<std::string> can_seqs;
	std::unordered_map<uint64_t, std::unordered_map<uint64_t, int>> edges;

	while(kseq_read(seq) >= 0){
		string name=seq->name.s;
		string sequence=seq->seq.s;
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		seq_len=seq->seq.l;

		string can_seq;
		uint64_t lastHash=0;
        std::unordered_set<uint64_t> seen_hashs;
		
		for (i = l = 0, x[0] = x[1] = 0; i < seq_len; ++i) {
			int c = seq_nt4_table[(uint8_t)sequence[i]];
			if (c < 4) {
				x[0] = (x[0] << 2 | c) & mask;
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= k) {
					string kmer=sequence.substr(i-k+1,k);
					uint64_t y = x[0] < x[1]? x[0] : x[1];
					uint64_t hash_value = hash64(y, mask);
					int Count=IntKmer2Value(h,hash_value);
					bool flag_kmer=checkKmer(kmer, Count, minCount);
                    bool flag_hash=seen_hashs.find(hash_value) == seen_hashs.end();
					
					if (flag_kmer && flag_hash){
						if (can_seq.length()==0){
							can_seq=kmer;
						}else{
							can_seq+=kmer.substr(k-1);
							edges[hash_value][lastHash]++;
						}
                        seen_hashs.insert(hash_value);
						lastHash=hash_value;
					}else{
						if (can_seq.length()>=(k+1)){
							can_seqs.push_back(can_seq);
						}
						can_seq="";
                        seen_hashs.clear();
					}
				}
			} 
			else {
				l = 0, x[0] = x[1] = 0;
			}
		}

        if (can_seq.length()>=(k+1)){
            can_seqs.push_back(can_seq);
        }
	}

	kseq_destroy(seq);
  	gzclose(fp);

	for (const auto& outerPair : edges) {
        std::vector<std::pair<uint64_t, int>> innerPairs(outerPair.second.begin(), outerPair.second.end());
        
        std::sort(innerPairs.begin(), innerPairs.end(),
                  [](const std::pair<uint64_t, int>& a, const std::pair<uint64_t, int>& b) {
                      return a.second > b.second;
                  });
        if (innerPairs[0].second<minCount){
            edges[outerPair.first].clear();
        }else{
            if (innerPairs.size()>=2){
                float ratio=static_cast<float>(innerPairs[0].second)/innerPairs[1].second;
                if (ratio<=2){
                    edges[outerPair.first].clear();
                } else {
                    auto& innerMap = edges[outerPair.first];
                    for (size_t i=1; i<innerPairs.size(); i++){
                        innerMap.erase(innerPairs[i].first);
                    }
                }
            }
        }
    }

    for (const auto& source : edges) {
        for (const auto& target : source.second) {
            std::cout << source.first << ":" << target.first << " : " << target.second<< std::endl;
        }
    }
    
    string best_adapter;
	int min_adapter_depth=0;
    int min_adapter_len=(P2In->AdapterLen);

    for (const auto& sequence : can_seqs) {
        string adapter;
		int adapter_depth=0;
        uint64_t lastHash;
        int seq_len=sequence.length();

        for (i = l = 0, x[0] = x[1] = 0; i < seq_len; ++i) {
			int c = seq_nt4_table[(uint8_t)sequence[i]];
            x[0] = (x[0] << 2 | c) & mask;                 
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
            if (++l >= k) {
                string kmer=sequence.substr(i-k+1,k);
                uint64_t y = x[0] < x[1]? x[0] : x[1];
                uint64_t hash_value = hash64(y, mask);

                if (edges.count(hash_value) > 0){
                    if (adapter.length()==0){
						adapter=kmer;
                    }else{
                        if(edges[hash_value].count(lastHash)>0){
                            adapter+=kmer.substr(k-1);
                            int depth=edges[hash_value][lastHash];
                            adapter_depth+=depth;

                        }else{
                            if (adapter.length()>=min_adapter_len){
                                if (adapter_depth>min_adapter_depth){
                                    best_adapter=adapter;
                                    min_adapter_depth=adapter_depth;
                                }
                            }
                            adapter=kmer;
                            adapter_depth=0;
                        }

                    }
                    lastHash=hash_value;
                }else{
                    if (adapter.length()>=min_adapter_len){
                        if (adapter_depth>min_adapter_depth){
                            best_adapter=adapter;
                            min_adapter_depth=adapter_depth;
                            cout <<best_adapter<<" "<<min_adapter_depth<<endl;
                        }
                    }
                    adapter="";
                    adapter_depth=0;
                }
            }
        }

        if (adapter.length()>=min_adapter_len){
            if (adapter_depth>min_adapter_depth){
                best_adapter=adapter;
                min_adapter_depth=adapter_depth;
                cout <<best_adapter<<" "<<min_adapter_depth<<endl;
            }
        }
    }

	mean_depth=int(min_adapter_depth/(best_adapter.length()));

	return best_adapter;
}

string Get_filter_parameter(Para_A24 *P2In, string &libname){

	string InPath=(P2In->InFile);

	gzFile fp;
  	kseq_t *seq;

	fp = gzopen(InPath.c_str(), "r");
  	seq = kseq_init(fp);

	int seqNum=0;
	int BCNum=P2In->BCNum;
	int ADNum=(P2In->ReadNumber);
	int maxSeq=BCNum;
	if (maxSeq<ADNum){
		maxSeq=ADNum;
	}
	int qtNum=5000;
	int minQ=50000;
	int maxQ=0;
	uint64_t sumQ = 0;
	uint64_t sumLen = 0;

	std::map<std::string, std::string> headSeq;
	std::map<std::string, std::string> tailSeq;

	int BCLen=P2In->BCLen;
	std::vector<int> headA(BCLen, 0);
	std::vector<int> headT(BCLen, 0);
	std::vector<int> tailA(BCLen, 0);
	std::vector<int> tailT(BCLen, 0);

	while((kseq_read(seq)) >= 0){
		if (seqNum>=maxSeq){
			break;
		}

		uint64_t length=seq->seq.l;
		if (length<500){
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

		string head_seq=sequence.substr(0, P2In->BCLen);
		string tail_seq=sequence.substr(length-(P2In->BCLen));
		if (seqNum<BCNum){
			for (int x=0; x<P2In->BCLen; x++){
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
			string head_ADSeq;
			string tail_ADSeq;
			if (P2In->BCLen!=P2In->EndLen){
				head_ADSeq=sequence.substr(0, P2In->EndLen);
				tail_ADSeq=sequence.substr(length-(P2In->EndLen));
			}else{
				head_ADSeq=head_seq;
				tail_ADSeq=tail_seq;
			}
			headSeq[name] = head_ADSeq;
			tailSeq[name] = tail_ADSeq;
		}
	}

	kseq_destroy(seq);
  	gzclose(fp);

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
	cout << "INFO: trim front length: "<<P2In->HeadCrop<<endl;
	cout << "INFO: trim tail length: "<<P2In->TailCrop<<endl;

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
		}else if(meanQ>=15){
			P2In->AverQ=10;
		}else if(meanQ>0) {
			P2In->AverQ=7;
		}else{
			P2In->AverQ=0;
		}
	}

	if ((P2In->Infq)==1){
		cout<<"INFO: min mean base quality was set to: "<<(P2In->AverQ)<<endl;
		P2In->AverQ=(P2In->AverQ)+qType;
	}

	//search adapter
	if (!(P2In->AdapterFile).empty() && access((P2In->AdapterFile).c_str(), 0) == 0){
		string adapter_name=P2In->AdapterFile;
		return adapter_name;
	}

	string prefix=GetFilePreifx(InPath);
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
	int depth_5p;
	int depth_3p;
	int len_5p=0;
	int len_3p=0;
	cout << "INFO: searching 5' adapter..."<<endl;
	adapter_5p=adapterSearch(P2In, libname, out_front, depth_5p);
	len_5p=adapter_5p.length();
	cout << "INFO: searching 3' adapter..."<<endl;
	adapter_3p=adapterSearch(P2In, libname, out_tail, depth_3p);
	len_3p=adapter_3p.length();

	if(len_5p==0 && len_3p==0){
		cout << "Warnings: the adapter could not be found within the general adapter library.";
		int p=10;
		uint64_t block_size = 10000000;

		cout << "INFO: assembly 5' adapter..."<<endl;
		kc_c4x_t *fh;
		fh = count_file(out_front.c_str(), P2In->Kmer, p, block_size, n_thread);
		adapter_5p=AssemblyAdapter(P2In, fh, out_front, depth_5p);
		len_5p=adapter_5p.length();
		for (int i = 0; i < 1<<p; ++i) kc_c4_destroy(fh->h[i]);
		free(fh->h); free(fh);

		cout << "INFO: assembly 3' adapter..."<<endl;
		kc_c4x_t *th;
		th = count_file(out_tail.c_str(), P2In->Kmer, p, block_size, n_thread);
		adapter_3p=AssemblyAdapter(P2In, th, out_tail, depth_3p);
		len_3p=adapter_3p.length();
		for (int i = 0; i < 1<<p; ++i) kc_c4_destroy(th->h[i]);
		free(th->h); free(th);
	}
	// write adapter
	string adapter_out;
	if(len_5p != 0 || len_3p != 0){
		adapter_out=prefix+".adapter.fa";
		ofstream OUTHA(adapter_out);
		if (len_5p >0 && len_3p >0){
			string rev_3p=rev_comp_seq(adapter_3p);
			if ((adapter_5p == adapter_3p) || (adapter_5p == rev_3p)){
				OUTHA <<">adapter_5-3p"<<"\n"<<adapter_5p<<endl;
			}else{
				OUTHA <<">adapter_5p"<<"\n"<<adapter_5p<<endl;
				OUTHA <<">adapter_3p"<<"\n"<<adapter_3p<<endl;
			}
		}else if(len_5p >0){
			OUTHA <<">adapter_5p"<<"\n"<<adapter_5p<<endl;
		}else if(len_3p >0){
			OUTHA <<">adapter_3p"<<"\n"<<adapter_3p<<endl;
		}
		OUTHA.close();
	}
	
	int minLen=len_5p;
	if (minLen>len_3p){
		minLen=len_3p;
	}
	P2In->MidLen=int(minLen*(P2In->Similarity));
	cout << "INFO: min match length in middle was set to: "<<P2In->MidLen<<endl;

	remove(out_front.c_str());
	remove(out_tail.c_str());
	return adapter_out;
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

void Split_low_qual_region_to_fastq(Para_A24 *P2In, string &raw_id, string &raw_seq, 
							string &raw_qual, int &seqLen, int &pass_num, string &OUT_DATA,
							std::array<uint64_t, 14>& DropInfo){
	string out_seq;
	string out_qual;
	string segment_seq;
	string segment_qual;
	double avgQuality;
	string seq_id=raw_id;
	int window=P2In->Window;
	int start_idx = 0;
	int end_idx = window;
	int lq_len=0;
	
	while (end_idx < seqLen) {
		segment_seq = raw_seq.substr(start_idx, window);
		segment_qual = raw_qual.substr(start_idx, window);
		avgQuality = CalcAvgQuality(segment_qual);
		if (avgQuality >= (P2In->AverQ)) {
			out_seq+=segment_seq;
			out_qual+=segment_qual;
		}else{
			lq_len+=window;
			if (out_seq.length()>=(P2In->MinLength)){
				pass_num++;
				DropInfo[2]++;
				DropInfo[3]+=out_seq.length();
				if (pass_num>=2){
					seq_id=raw_id+":"+std::to_string(pass_num);
				}
				OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + out_seq + 
							"\n+\n" + out_qual + "\n";
			}else{
				if (out_seq.length()>0){
					DropInfo[12]++;
					DropInfo[13]+=out_seq.length();
				}
			}
			out_seq="";
			out_qual="";
		}
		start_idx+=window;
		end_idx=start_idx+window;
	}
	//Processing the final segment
	segment_seq = raw_seq.substr(start_idx);
	segment_qual = raw_qual.substr(start_idx);
	avgQuality = CalcAvgQuality(segment_qual);
	if (avgQuality >= (P2In->AverQ)){
		out_seq+=segment_seq;
		out_qual+=segment_qual;
	}else{
		lq_len+=(seqLen-start_idx);
	}

	if (out_seq.length()>=(P2In->MinLength)){
		pass_num++;
		DropInfo[2]++;
		DropInfo[3]+=out_seq.length();
		if (pass_num>=2){
			seq_id=raw_id+":"+std::to_string(pass_num);
		}
		OUT_DATA=OUT_DATA+"@"+seq_id+"\n" + out_seq + 
					"\n+\n" + out_qual + "\n";
	}else{
		if (out_seq.length()>0){
			DropInfo[12]++;
			DropInfo[13]+=out_seq.length();
		}
	}

	if (lq_len>0){
		DropInfo[10]++;
		DropInfo[11]+=lq_len;
	}
}

void Split_low_qual_region_to_fasta(Para_A24 *P2In, string &raw_id, string &raw_seq, 
							string &raw_qual, int &seqLen, int &pass_num, string &OUT_DATA,
							std::array<uint64_t, 14>& DropInfo){
	string out_seq;
	string segment_seq;
	string segment_qual;
	double avgQuality;
	string seq_id=raw_id;
	int window=P2In->Window;
	int start_idx = 0;
	int end_idx = window;
	int lq_len=0;
	
	while (end_idx < seqLen) {
		segment_seq = raw_seq.substr(start_idx, window);
		segment_qual = raw_qual.substr(start_idx, window);
		avgQuality = CalcAvgQuality(segment_qual);
		if (avgQuality >= (P2In->AverQ)) {
			out_seq+=segment_seq;
		}else{
			lq_len+=window;
			if (out_seq.length()>=(P2In->MinLength)){
				pass_num++;
				DropInfo[2]++;
				DropInfo[3]+=out_seq.length();
				if (pass_num>=2){
					seq_id=raw_id+":"+std::to_string(pass_num);
				}
				OUT_DATA=OUT_DATA+">"+seq_id+"\n" + out_seq + "\n";
			}else{
				if (out_seq.length()>0){
					DropInfo[12]++;
					DropInfo[13]+=out_seq.length();
				}
			}
			out_seq="";
		}
		start_idx+=window;
		end_idx=start_idx+window;
	}
	//Processing the final segment
	segment_seq = raw_seq.substr(start_idx);
	segment_qual = raw_qual.substr(start_idx);
	avgQuality = CalcAvgQuality(segment_qual);
	if (avgQuality >= (P2In->AverQ)){
		out_seq+=segment_seq;
	}else{
		lq_len+=(seqLen-start_idx);
	}

	if (out_seq.length()>=(P2In->MinLength)){
		pass_num++;
		DropInfo[2]++;
		DropInfo[3]+=out_seq.length();
		if (pass_num>=2){
			seq_id=raw_id+":"+std::to_string(pass_num);
		}
		OUT_DATA=OUT_DATA+">"+seq_id+"\n" + out_seq + "\n";
	} else{
		if (out_seq.length()>0){
			DropInfo[12]++;
			DropInfo[13]+=out_seq.length();
		}
	}

	if (lq_len>0){
		DropInfo[10]++;
		DropInfo[11]+=lq_len;
	}
}

void adapterMap(Para_A24 * P2In, string &raw_seq, int &raw_len,
				const mm_idx_t* mi, mm_mapopt_t &mopt,
				const mm_idx_t* emi, mm_mapopt_t &emopt, mm_tbuf_t *tbuf,
				std::vector<std::vector<int>> &keepRegions,
				std::array<uint64_t, 14>& DropInfo){

	std::vector<std::vector<int>> adapterRegions;
	int minLen=UINT64_MAX;
	
	//for middle
	mm_reg1_t *reg;
	int j, n_reg;
	reg = mm_map(mi, raw_len, raw_seq.c_str(), &n_reg, tbuf, &mopt, 0);
	for (j = 0; j < n_reg; ++j) {
		mm_reg1_t *r = &reg[j];
		int ql=raw_len;
		int qs=r->qs;
		int qe=r->qe;
		int tl=mi->seq[r->rid].len;
		int ts=r->rs;
		int te=r->re;

		int mlen=r->mlen;
		if(mlen<P2In->MidLen){
			continue;
		}

		int blen=r->blen;
		if (blen<tl){
			blen=tl;
		}

		if (mlen<(blen*(P2In->Similarity))){
			continue;
		}
		if (tl<minLen){
			minLen=tl;
		}
		adapterRegions.push_back({qs, qe});
	}
	free(reg);
	//
	int endLen=minLen;
	if (endLen<50){
		endLen=50;
	}

	//for head end
	string head_seq=raw_seq.substr(0,endLen);
	int head_len=endLen;
	mm_reg1_t *hreg;
	int hj, hn_reg;
	hreg = mm_map(emi, head_len, head_seq.c_str(), &hn_reg, tbuf, &emopt, 0);
	
	for (hj = 0; hj < hn_reg; ++hj) {
		mm_reg1_t *r = &hreg[hj];
		int ql=head_len;
		int qs=r->qs;
		int qe=r->qe;
		int tl=emi->seq[r->rid].len;
		int ts=r->rs;
		int te=r->re;
		int mlen=r->mlen;

		if(mlen<P2In->MatchLen){
			continue;
		}

		int blen=qe;
		if (mlen<(blen*(P2In->Similarity))){
			continue;
		}
		adapterRegions.push_back({qs, qe});
	}
	free(hreg);

	//for tail end
	string tail_seq=raw_seq.substr(raw_len-endLen,endLen);
	int tail_len=endLen;
	mm_reg1_t *treg;
	int tj, tn_reg;
	treg = mm_map(emi, tail_len, tail_seq.c_str(), &tn_reg, tbuf, &emopt, 0);
	for (tj = 0; tj < tn_reg; ++tj) {
		mm_reg1_t *r = &treg[tj];
		int ql=tail_len;
		int qs=r->qs;
		int qe=r->qe;
		int tl=emi->seq[r->rid].len;
		int ts=r->rs;
		int te=r->re;
		int mlen=r->mlen;

		if(mlen<P2In->MatchLen){
			continue;
		}
		int blen=ql-qs;
		if (mlen<(blen*(P2In->Similarity))){
			continue;
		}
		adapterRegions.push_back({qs+raw_len-endLen, qe+raw_len-endLen});
	}
	free(treg);

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

	if (mergedRegions.size()>=1){
		DropInfo[8]++;
		int currentStart = 0;
		for (const auto& region : mergedRegions) {
			DropInfo[9]+=(region[1]-region[0]);
			if (region[0] > currentStart) {
				keepRegions.push_back({currentStart, region[0] - currentStart});
			}
			currentStart = region[1];
		}

		if (currentStart < raw_len) {
			keepRegions.push_back({currentStart, static_cast<int>(raw_len) - currentStart});
		}
	}
}

void Filter_fastq_reads_adapter(Para_A24 * P2In, string &OUT_DATA,
						vector<string>& ID, vector <string> &SEQ, vector <string> &QUAL,
						const mm_idx_t* mi, mm_mapopt_t &mopt, 
						const mm_idx_t* emi, mm_mapopt_t &emopt,mm_tbuf_t *tbuf,
						std::array<uint64_t, 14>& DropInfo){

	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);
	
	for (int i = 0; i < SEQ.size(); i++) {
		int seq_len = SEQ[i].length();
		int qual_len = QUAL[i].length();
		DropInfo[0]++;
		DropInfo[1]+=seq_len;

		if ((totalCrop >= seq_len) || (seq_len != qual_len) || (seq_len >(P2In->MaxLength))){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
			continue;
		}

		string raw_seq = SEQ[i].substr(headCrop, seq_len - totalCrop);
		string raw_qual = QUAL[i].substr(headCrop, qual_len - totalCrop);
		int raw_len = raw_seq.length();

		if (raw_len < (P2In->MinLength)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
			continue;
		}

		if (totalCrop>0){
			DropInfo[6]++;
			DropInfo[7]+=totalCrop;
		}

		std::vector<std::vector<int>> keepRegions;
		adapterMap(P2In, raw_seq, raw_len, mi, mopt, emi, emopt, tbuf, keepRegions,DropInfo);

		string raw_id=ID[i];
		int pass_num=0;
		if ((P2In->Outfq)==1){
			if (keepRegions.size()==0){
				Split_low_qual_region_to_fastq(P2In, raw_id, raw_seq, raw_qual, 
											  raw_len, pass_num, OUT_DATA,DropInfo);
			}else{
				for (const auto& region : keepRegions) {
					int start = region[0];
					int seqLen = region[1];
					string seq=raw_seq.substr(start,seqLen);
					string qual=raw_qual.substr(start,seqLen);

					Split_low_qual_region_to_fastq(P2In, raw_id, seq, qual, 
												seqLen, pass_num, OUT_DATA,DropInfo);
				}
			}
		}else{
			if (keepRegions.size()==0){
				Split_low_qual_region_to_fasta(P2In, raw_id, raw_seq, raw_qual, 
									raw_len, pass_num, OUT_DATA,DropInfo);
			}else{
				for (const auto& region : keepRegions) {
					int start = region[0];
					int seqLen = region[1];
					string seq=raw_seq.substr(start,seqLen);
					string qual=raw_qual.substr(start,seqLen);

					Split_low_qual_region_to_fasta(P2In, raw_id, seq, qual, 
												seqLen, pass_num, OUT_DATA,DropInfo);
				}
			}
		}
	}
}

void Filter_fastq_reads(Para_A24 * P2In, string &OUT_DATA, vector<string>& ID, 
						vector <string> &SEQ, vector <string> &QUAL,
						std::array<uint64_t, 14>& DropInfo){
	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);

	for (int i = 0; i < SEQ.size(); i++) {
		int seq_len = SEQ[i].length();
		int qual_len = QUAL[i].length();
		DropInfo[0]++;
		DropInfo[1]+=seq_len;
		
		if ((totalCrop >= seq_len) || (seq_len != qual_len) || (seq_len >(P2In->MaxLength))) {
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
			continue;
		}

		string raw_seq = SEQ[i].substr(headCrop, seq_len - totalCrop);
		string raw_qual = QUAL[i].substr(headCrop, qual_len - totalCrop);
		int raw_len = raw_seq.length();
		
		if (raw_len < (P2In->MinLength)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
			continue;
		}

		if (totalCrop>0){
			DropInfo[6]++;
			DropInfo[7]+=totalCrop;
		}

		string raw_id=ID[i];
		int pass_num=0;
		if ((P2In->Outfq)==1){
			Split_low_qual_region_to_fastq(P2In, raw_id, raw_seq, raw_qual, 
											raw_len, pass_num, OUT_DATA,DropInfo);
		}else{
			Split_low_qual_region_to_fasta(P2In, raw_id, raw_seq, raw_qual, 
											raw_len, pass_num, OUT_DATA,DropInfo);
		}
	}    
}

void Filter_fasta_reads_adapter(Para_A24 * P2In,string &OUT_DATA, 
						vector<string>& ID, vector <string> & SEQ,
                        const mm_idx_t* mi, mm_mapopt_t &mopt, 
						const mm_idx_t* emi, mm_mapopt_t &emopt, mm_tbuf_t *tbuf,
						std::array<uint64_t, 14>& DropInfo) {
	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);
	
	for (int i = 0; i < SEQ.size(); i++) {
		int seq_len = SEQ[i].length();
		DropInfo[0]++;
		DropInfo[1]+=seq_len;

        if (totalCrop >= seq_len || seq_len >(P2In->MaxLength)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
            continue;
        }

        string raw_seq = SEQ[i].substr(headCrop, seq_len - totalCrop);
        int raw_len = raw_seq.length();

        if (raw_len < (P2In->MinLength)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
            continue;
        }

		if (totalCrop>0){
			DropInfo[6]++;
			DropInfo[7]+=totalCrop;
		}

        std::vector<std::vector<int>> keepRegions;
		adapterMap(P2In, raw_seq, raw_len, mi, mopt, emi, emopt, tbuf, keepRegions,DropInfo);
        
        string seq_id=ID[i];
        int pass_num=0;
        for (const auto& region : keepRegions) {
            int start = region[0];
            int seqLen = region[1];
            if (seqLen < (P2In->MinLength)){
				DropInfo[12]++;
				DropInfo[13]+=seqLen;
				continue;
			}

            pass_num++;
			DropInfo[2]++;
			DropInfo[3]+=seqLen;

            if (pass_num>=2){
                seq_id=ID[i]+":"+std::to_string(pass_num);
            }
            string seq=raw_seq.substr(start,seqLen);
            OUT_DATA=OUT_DATA+">"+seq_id+"\n" + seq + "\n";
        }
	}
}

void Filter_fasta_reads(Para_A24 * P2In,string &OUT_DATA, vector<string> &ID, 
						vector <string> &SEQ, std::array<uint64_t, 14>& DropInfo) {
	int headCrop = P2In->HeadCrop;
	int totalCrop = (P2In->HeadCrop) + (P2In->TailCrop);

	for (int i = 0; i < SEQ.size(); i++) {
		int seq_len = SEQ[i].length();
		DropInfo[0]++;
		DropInfo[1]+=seq_len;

		if (totalCrop >= seq_len || seq_len >(P2In->MaxLength)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
            continue;
        }

		string seq = SEQ[i].substr(headCrop, seq_len - totalCrop);
		int seqLen = seq.length();
		if (seqLen < (P2In->MinLength)){
			DropInfo[4]++;
			DropInfo[5]+=seq_len;
			continue;
		}

		if (totalCrop>0){
			DropInfo[6]++;
			DropInfo[7]+=totalCrop;
		}

		DropInfo[2]++;
		DropInfo[3]+=seqLen;

		OUT_DATA=OUT_DATA+">"+ID[i]+"\n" + seq + "\n";
	}
}

void Filter_reads_adapter(Para_A24 * P2In, string &OUT_DATA, 
				vector<string> &ID, vector <string> &SEQ, vector <string> &QUAL,
				const mm_idx_t* mi, mm_mapopt_t &mopt, 
				const mm_idx_t* emi, mm_mapopt_t &emopt, mm_tbuf_t *tbuf,
				std::array<uint64_t, 14>& DropInfo){
	OUT_DATA="";
	if(QUAL[0].length()<=2){
		Filter_fasta_reads_adapter(P2In, OUT_DATA, ID, SEQ, mi, mopt,emi, emopt, tbuf,
								   DropInfo);
	} else {
		Filter_fastq_reads_adapter(P2In, OUT_DATA, ID, SEQ, QUAL, mi, mopt, emi, emopt,tbuf,
								   DropInfo);
	}
}

void Filter_reads(Para_A24 * P2In, string &OUT_DATA, vector<string> &ID, 
				vector <string> &SEQ, vector <string> &QUAL, std::array<uint64_t, 14>& DropInfo){
	OUT_DATA="";
	if(QUAL[0].length()<=2){
		Filter_fasta_reads(P2In, OUT_DATA, ID, SEQ, DropInfo);
	} else {
		Filter_fastq_reads(P2In, OUT_DATA, ID, SEQ, QUAL, DropInfo);
	}
}

void Filter_reads_adapter_gz(Para_A24 *P2In, vector<string> &ID, vector <string> &SEQ, 
					vector <string> &QUAL, const mm_idx_t* mi, mm_mapopt_t &mopt,
					const mm_idx_t* emi, mm_mapopt_t &emopt, mm_tbuf_t *tbuf,
					std::array<uint64_t, 14>& DropInfo,
					uint8_t ** ComData, size_t & ComSize,
					size_t &ComBuff, int &tid){
	string OUT_DATA="";
	OUT_DATA.reserve(OUTPUT_BUFFER_SIZE);
	if(QUAL[0].length()<=2){
		Filter_fasta_reads_adapter(P2In, OUT_DATA, ID, SEQ, mi, mopt,emi, emopt, tbuf, 
									DropInfo);

	} else {
		Filter_fastq_reads_adapter(P2In, OUT_DATA, ID, SEQ, QUAL, mi, mopt,emi, emopt,tbuf,
									DropInfo);
		
	}
	
	if (!(OUT_DATA.empty())) {
		DeflateCompress  GZData;
		size_t inputSize  = OUT_DATA.length();
		GZData.compressData(OUT_DATA.c_str(), inputSize, ComData, 
							ComSize,ComBuff,tid);
	}
}

void Filter_reads_gz(Para_A24 * P2In, vector<string> &ID, vector <string> &SEQ, 
					vector <string> &QUAL, std::array<uint64_t, 14>& DropInfo,
					uint8_t ** ComData, size_t & ComSize,
					size_t & ComBuff, int &tid){
	string OUT_DATA="";
	OUT_DATA.reserve(OUTPUT_BUFFER_SIZE);
	if(QUAL[0].length()<=2){
		Filter_fasta_reads(P2In, OUT_DATA, ID, SEQ, DropInfo);

	} else {
		Filter_fastq_reads(P2In, OUT_DATA, ID, SEQ, QUAL, DropInfo);
		
	}
	
	if (!(OUT_DATA.empty())) {
		DeflateCompress  GZData;
		size_t inputSize  = OUT_DATA.length();
		GZData.compressData(OUT_DATA.c_str(), inputSize, ComData, 
							ComSize,ComBuff,tid);
	}
}

int Run_seq_filter_adapter(Para_A24 * P2In,string &adapter_name, 
						   mm_idxopt_t &iopt, mm_mapopt_t &mopt,
						   mm_idxopt_t &eiopt, mm_mapopt_t &emopt) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	std::vector<std::vector<std::string>> ID(n_thread);
	std::vector<std::vector<std::string>> SEQ(n_thread);
	std::vector<std::vector<std::string>> QUAL(n_thread);

	std::vector<std::array<uint64_t, 14>> DropInfo(n_thread, std::array<uint64_t, 14>{});
	//raw_reads raw_bases
	//clean_reads clean_bases
	//shortDrop_reads shortDrop_bases
	//endDrop_reads  endDrop_bases
	//adapterDrop_reads adapterDrop_bases
	//LQDrop_reads LQDrop_bases
	//OutDrop_reads OutDrop_bases

	std::vector<std::thread> threads;

	gzFile fp = gzopen((P2In->InFile).c_str(), "r");
    kseq_t *seq = kseq_init(fp);

	string outname=P2In->OutFile;
	int total_length = 0;
	//
	int index_thread=1;
	//
	mm_tbuf_t *tbufs[n_thread];
	for (int i = 0; i < n_thread; i++) {
		tbufs[i] = mm_tbuf_init();
	}

	if (P2In->OUTGZ){
		DeflateOgzstream OUTHGZ(outname.c_str());
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
		mm_idx_reader_t *r = mm_idx_reader_open(adapter_name.c_str(), &iopt, 0);
		mm_idx_t *mi;

		mm_idx_reader_t *er = mm_idx_reader_open(adapter_name.c_str(), &eiopt, 0);
		mm_idx_t *emi;

		while ((mi = mm_idx_reader_read(r, index_thread)) != 0 && 
			   (emi = mm_idx_reader_read(er, index_thread)) != 0) {

			mm_mapopt_update(&mopt, mi);
			mm_mapopt_update(&emopt, emi);

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
										ref(SEQ[tid]),ref(QUAL[tid]),ref(mi),ref(mopt),ref(emi),ref(emopt),
										ref(tbufs[tid]),ref(DropInfo[tid]),
										ComData, ref(ComSize[tid]),
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
								OUTHGZ.writeGZIO(ComData[i], ComSize[i] );
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
								ref(SEQ[tid]),ref(QUAL[tid]),ref(mi), ref(mopt),ref(emi),ref(emopt),
								ref(tbufs[tid]),ref(DropInfo[tid]),
								ComData, ref(ComSize[tid]),
								ref(ComBuff[tid]), ref(ArryThread[tid])));
				total_length = 0;
			}
			
			for (auto &t : threads) {
				t.join();
			}

			threads.clear();

			for (int i = 0; i <= tid; i++) {
				if (ComSize[i]>0) {
					OUTHGZ.writeGZIO(ComData[i], ComSize[i] );
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

        	mm_idx_destroy(mi);
		}
		mm_idx_reader_close(r);

	} else {
		ofstream OUTH;
		OUTH.open(outname.c_str());
		vector <string> OUT_DATA;
		OUT_DATA.resize(n_thread);

		int tid=0;
		mm_idx_reader_t *r = mm_idx_reader_open(adapter_name.c_str(), &iopt, 0);
		mm_idx_t *mi;

		mm_idx_reader_t *er = mm_idx_reader_open(adapter_name.c_str(), &eiopt, 0);
		mm_idx_t *emi;

		while ((mi = mm_idx_reader_read(r, index_thread)) != 0 && 
			   (emi = mm_idx_reader_read(er, index_thread)) != 0) {

			mm_mapopt_update(&mopt, mi);
			mm_mapopt_update(&emopt, emi);

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
							ref(ID[tid]),ref(SEQ[tid]),ref(QUAL[tid]),
							ref(mi), ref(mopt),ref(emi),ref(emopt),ref(tbufs[tid]),
							ref(DropInfo[tid])));
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
						ref(ID[tid]),ref(SEQ[tid]),ref(QUAL[tid]),
						ref(mi), ref(mopt),ref(emi),ref(emopt),ref(tbufs[tid]),
						ref(DropInfo[tid])));
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
        	mm_idx_destroy(mi);
		}
		mm_idx_reader_close(r);
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

	for (int i = 0; i < n_thread; i++) {
		mm_tbuf_destroy(tbufs[i]);
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
	}

	cout << raw_reads<<" reads with "<<raw_bases<<" bases were input"<<endl;
	cout << shortDrop_reads <<" reads with "<<shortDrop_bases<<" bases were dropped before filter"<<endl;
	cout << endDrop_reads <<" reads were trimmed "<<endDrop_bases<<" bases in end"<<endl;
	cout << adapterDrop_reads <<" reads were trimmed " <<adapterDrop_bases<<" bases with adapter"<<endl;
	cout << LQDrop_reads <<" reads were trimmed "<<LQDrop_bases<<" bases with low quality regions"<<endl;
	cout << outDrop_reads<<" reads with "<<outDrop_bases<<" bases were dropped before output"<<endl;
	cout << clean_reads <<" reads with "<<clean_bases<<" bases were output"<<endl;

	return 0;
}

int Run_seq_filter(Para_A24 * P2In) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	std::vector<std::vector<std::string>> ID(n_thread);
	std::vector<std::vector<std::string>> SEQ(n_thread);
	std::vector<std::vector<std::string>> QUAL(n_thread);

	std::vector<std::array<uint64_t, 14>> DropInfo(n_thread, std::array<uint64_t, 14>{});
	//raw_reads raw_bases
	//clean_reads clean_bases
	//shortDrop_reads shortDrop_bases
	//endDrop_reads  endDrop_bases
	//adapterDrop_reads adapterDrop_bases
	//LQDrop_reads LQDrop_bases
	//OutDrop_reads OutDrop_bases

	std::vector<std::thread> threads;

	gzFile fp = gzopen((P2In->InFile).c_str(), "r");
    kseq_t *seq = kseq_init(fp);

	string outname=P2In->OutFile;
	//

	int total_length = 0;
	if (P2In->OUTGZ){
		DeflateOgzstream OUTHGZ(outname.c_str());
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
				threads.push_back(thread(Filter_reads_gz,P2In, ref(ID[tid]), 
									ref(SEQ[tid]),ref(QUAL[tid]),ref(DropInfo[tid]),
									ComData, ref(ComSize[tid]),
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
							OUTHGZ.writeGZIO(ComData[i], ComSize[i] );
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
			threads.push_back(thread(Filter_reads_gz,P2In, ref(ID[tid]), 
							ref(SEQ[tid]),ref(QUAL[tid]),ref(DropInfo[tid]),
							ComData, ref(ComSize[tid]),
							ref(ComBuff[tid]), ref(ArryThread[tid])));
			total_length = 0;
		}
		
		for (auto &t : threads) {
			t.join();
		}

		threads.clear();

		for (int i = 0; i <= tid; i++) {
			if (ComSize[i]>0) {
				OUTHGZ.writeGZIO(ComData[i], ComSize[i] );
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
                
				threads.push_back(thread(Filter_reads,P2In, ref(OUT_DATA[tid]), 
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
			threads.push_back(thread(Filter_reads,P2In, ref(OUT_DATA[tid]), 
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
	}
	cout << raw_reads<<" reads with "<<raw_bases<<" bases were input"<<endl;
	cout << shortDrop_reads <<" reads with "<<shortDrop_bases<<" bases were dropped before filter"<<endl;
	cout << endDrop_reads <<" reads were trimmed "<<endDrop_bases<<" bases in end"<<endl;
	cout << adapterDrop_reads <<" reads were trimmed " <<adapterDrop_bases<<" bases with adapter"<<endl;
	cout << LQDrop_reads <<" reads were trimmed "<<LQDrop_bases<<" bases with low quality regions"<<endl;
	cout << outDrop_reads<<" reads with "<<outDrop_bases<<" bases were dropped before output"<<endl;
	cout << clean_reads <<" reads with "<<clean_bases<<" bases were output"<<endl;

	return 0;
}

int Run_TGSFilter(Para_A24 * P2In, string &libname) {
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
	string adapter_name=Get_filter_parameter(P2In,libname);

	if (P2In->ONLYAD){
		return 0;
	}

	//
	int mk=P2In->MatchLen;
	if (mk % 2 != 1) {
        mk--;
    }

	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	int cnt=(P2In->MidLen)/mk;
	mm_verbose = 2;
	mm_set_opt(0, &iopt, &mopt);
	iopt.k=15;
	iopt.w=1;
	mopt.flag |= MM_F_CIGAR;
	mopt.min_cnt = cnt;
	mopt.min_chain_score =mk;
	mopt.min_dp_max=mk;
	//
	
	mm_idxopt_t eiopt;
	mm_mapopt_t emopt;

	mm_verbose = 2;
	mm_set_opt(0, &eiopt, &emopt);
	eiopt.k=mk;
	eiopt.w=1;
	emopt.flag |= MM_F_CIGAR;
	emopt.min_cnt = 1;
	emopt.min_chain_score =mk;
	emopt.min_dp_max=mk;

	if (adapter_name.empty()){
		Run_seq_filter(P2In);
	}else{
		Run_seq_filter_adapter(P2In,adapter_name, iopt, mopt, eiopt, emopt);
	}

	return 0;
}

std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string get_lib_path(const std::string& programPath) {
    std::string check_path;

    size_t found = programPath.find_last_of("/\\");
    if (found != std::string::npos) {
        check_path = programPath.substr(0, found) + "/tgs_adapter.fa";
        if (access(check_path.c_str(), F_OK) == 0) {
            return check_path;
        }
    }

    const char* pathEnv = std::getenv("PATH");
    if (pathEnv == nullptr) {
        return "";
    }

    std::vector<std::string> paths = split(pathEnv, ':');
    for (const auto& path : paths) {
        check_path = path + "/tgs_adapter.fa";
        if (access(check_path.c_str(), F_OK) == 0) {
            return check_path;
        }
    }

    return "";
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

	std::string programPath = argv[0];
	std::string libname=get_lib_path(programPath);
	//cout << "lib name: "<<libname<<endl;

	Run_TGSFilter(P2In,libname);

	delete P2In ;
	return 0;
}

