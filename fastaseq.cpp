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
	int header_offset;
	int seq_offset;
	int seq_length;
	int line_width;
};


string INDEX_FILE_EXT = ".idx";
string SEP = "\t";


vector<index_entry> createIndex(string filename) {
	ifstream fasta;
	fasta.open(filename.c_str());
	if(!fasta.is_open()) {
		cerr << "could not open file " << filename << endl;
	    exit(1);
	}
	ofstream idx;
	idx.open((filename + INDEX_FILE_EXT).c_str());
	if(!idx.is_open()) {
		fasta.close();
		cerr << "could not create index file " << filename << endl;
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
				idx << e.title << SEP << e.header_offset << SEP << e.seq_offset << SEP << e.seq_length << SEP << e.line_width << endl;
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
		idx << e.title << SEP << e.header_offset << SEP << e.seq_offset << SEP << e.seq_length << SEP << e.line_width << endl;
		list.push_back(e);
	}
	
	idx.flush();
	idx.close();
	
	return list;
}


vector<index_entry> readIndex(string filename) {
	ifstream file;
	file.open((filename + INDEX_FILE_EXT).c_str());
	if(!file.is_open()) {
		return createIndex(filename);
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
	map<char,char> tr = {{'A','T'}, {'T','A'}, {'G','C'}, {'C','G'}, {'N','N'}};
	for(int i = seq.length()-1; i >= 0; i--) {
		rev << tr[seq[i]];
	}
	return rev.str();
}


string getRandomSequence(FILE* fasta, vector<index_entry> index, int length) {
	index_entry e = index[rand() % index.size()];
	int pos = -1;
	while(pos < 0 || pos > e.seq_length - length) {
		pos = rand() % e.seq_length;
	}
	
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


int main(const int argc, const char** argv) {
	vector<index_entry> index = readIndex(string(argv[1]));
	FILE* fasta = fopen(argv[1], "r");
	srand(time(0));
	cerr << getRandomSequence(fasta, index, 99) << endl;
	cerr << getRandomSequence(fasta, index, 99) << endl;
	cerr << getRandomSequence(fasta, index, 99) << endl;
	cerr << getRandomSequence(fasta, index, 99) << endl;
	fclose(fasta);
}


