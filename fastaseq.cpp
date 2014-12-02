#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;


struct index_entry {
	string title;
	unsigned int header_offset;
	unsigned int seq_offset;
	unsigned int seq_length;
	unsigned int line_width;
};


const string INDEX_FILE_EXT = ".idx";
const string SEP = "\t";
static map<char,char> COMPLEMENT = {{'A','T'}, {'T','A'}, {'G','C'}, {'C','G'}, {'N','N'}};


vector<index_entry> createIndex(string fastaname, string indexname) {
	ifstream fasta;
	fasta.open(fastaname.c_str());
	if(!fasta.is_open()) {
		cerr << "could not open file " << fastaname << endl;
	    exit(1);
	}
	ofstream idx;
	idx.open(indexname.c_str());
	if(!idx.is_open()) {
		fasta.close();
		cerr << "could not create index file " << indexname << endl;
	    exit(1);
	}
	string line;
	vector<index_entry> list;
	bool first_entry = true;
	index_entry e;
	int offset = 0;
	while(getline(fasta, line)) {
		if(line[0] == '>' || line[0] == '@') {
			if(!first_entry) {
				idx << e.title << SEP << e.header_offset << SEP << e.seq_offset;
				idx << SEP << e.seq_length << SEP << e.line_width << endl;
				list.push_back(e);
			} else {
				first_entry = false;
			}
			e = index_entry();
			e.title = line.substr(1,line.find_first_of(" \t\n")-1);
			e.header_offset = offset;
			e.seq_offset = offset + line.length() + 1;
		} else if(line.length() == 0 || line[0] == ';' || line[0] == '+' || line[0] == ' ') {
			if(e.seq_length == 0) e.seq_offset = offset + line.length() + 1;
		} else {
			e.seq_length += line.length();
			if(line.length() > e.line_width) e.line_width = line.length();
		}
		offset += line.length() + 1;
	}
	fasta.close();
	
	if(e.seq_offset > 0) {
		idx << e.title << SEP << e.header_offset << SEP << e.seq_offset; 
		idx << SEP << e.seq_length << SEP << e.line_width << endl;
		list.push_back(e);
	}
	
	idx.flush();
	idx.close();
	
	return list;
}


vector<index_entry> readIndex(string fastaname, string indexname) {
	ifstream file;
	file.open(fastaname.c_str());
	if(!file.good()) {
		return createIndex(fastaname, indexname);
	}
	vector<index_entry> index;
	string line;
	while(getline(file, line)) {
		istringstream s(line);
		index_entry e;
		s >> e.title;
		s >> e.header_offset;
		s >> e.seq_offset;
		s >> e.seq_length;
		s >> e.line_width;
		index.push_back(e);
	}
	file.close();
	return index;
}


string reverse_complement(string seq) {
	stringstream rev;
	for(int i = seq.length()-1; i >= 0; i--) {
		rev << COMPLEMENT[seq[i]];
	}
	return rev.str();
}


string getSubsequence(FILE* fasta, index_entry e, unsigned int from, unsigned int to) {
	unsigned int pos = from;
	unsigned int length = to-from;
	bool minus_strand = false;
	if(from > to) { pos = to; length = from-to; minus_strand = true; }
	int num_newlines = (e.line_width - (pos % e.line_width) + length) / e.line_width;
    int seqlen = length + num_newlines;
    char* seq = (char*) calloc (seqlen + 1, sizeof(char));
    int offset = e.seq_offset + pos + (pos / e.line_width);
    fseek(fasta, offset , SEEK_SET);
    fread(seq, sizeof(char), seqlen, fasta);
    seq[seqlen] = '\0';
    char* pbegin = seq;
    char* pend = seq + (seqlen/sizeof(char));
    pend = remove(pbegin, pend, '\n');
    pend = remove(pbegin, pend, '\0');
    string s = seq;
    free(seq);
    s.resize((length)/sizeof(char));
    if(minus_strand) s = reverse_complement(s);
    return s;
}


string getRandomSequence(FILE* fasta, vector<index_entry> index, unsigned int length) {
	index_entry e = index[rand() % index.size()];
	int pos = -1;
	while(pos < 0 || pos > (int)e.seq_length - (int)length) {
		pos = rand() % e.seq_length;
	}
	string s = getSubsequence(fasta, e, pos, pos+length);
	stringstream header;
	header << ">" << e.title;
	if(rand() % 2 == 1) {
    		s = reverse_complement(s);
		header << " complement " << (pos+length+1) << "-" << (pos+1) << endl;
	} else {
		header << " " << (pos+1) << "-" << (pos+length+1) << endl;
	}
	header << s;
	return header.str();
}


void print_help() {
	cout << "This is fastaseq 0.2";
}


int main(const int argc, const char** argv) {
	string fastaname; string indexname; int len_seq = 100; int batch = 1;
	for(int i=1; i<argc; i++) {
		string p = string(argv[i]);
		if(p.compare("-h") == 0 || p.compare("--help") == 0) { 
			print_help(); 
		} else if(p.compare("-f") == 0 || p.compare("--fasta") == 0) {
			fastaname = string(argv[++i]);
		} else if(p.compare("-i") == 0 || p.compare("--index") == 0) {
			indexname = string(argv[++i]);
		} else if(p.compare("-l") == 0 || p.compare("--length") == 0) {
			istringstream ss(argv[++i]);
			if(!(ss >> len_seq)) { cerr << "illegal length (-l) parameter: " << argv[i] << endl; return -1; }
		} else if(p.compare("-b") == 0 || p.compare("--batch") == 0) {
			istringstream ss(argv[++i]);
			if(!(ss >> batch)) { cerr << "illegal batch (-b) size: " << argv[i] << endl; return -1; }
		}
	}
	if(fastaname.empty()) { 
		cerr << "no fasta file parameter was given! use '-f' option" << endl; 
		print_help();
		return -1;
	}
	if(indexname.empty()) { indexname = fastaname + INDEX_FILE_EXT; createIndex(fastaname, indexname); }
	vector<index_entry> index = readIndex(fastaname, indexname);
	FILE* fasta = fopen(fastaname.c_str(), "r");
	srand(time(0));
	for(int i = 0; i<batch; i++) {
		cout << getRandomSequence(fasta, index, len_seq) << endl;
	}
	cerr << getRandomSequence(fasta, index, 99) << endl;
	fclose(fasta);
}


