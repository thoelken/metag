#include <string>
#include <iostream>
#include <array>
#include <vector>
#include <bitset>
#include <cmath>

//const double M_PI = 3.14159265358979323846;

using namespace std;

/// encode capital nucleotide character as bitset
bitset<4> encode_nt(char nt) {
	bitset<4> encoding;// = {0,0,0,0};
	switch(nt) {
		case 'A': encoding[0] = 1; break;
		case 'C': encoding[1] = 1; break;
		case 'G': encoding[2] = 1; break;
		case 'T': encoding[3] = 1; break;
		default: cerr << "nucleotide unaccepted in sequence " << nt << endl;
	}
	return encoding;
}

/// Discrete Cosine Transformation of binary signal. Returns only N coefficients (default=5)
vector<double> dct(vector<bool> x, int N = 5) {
	vector<double> output;
	if(N < 1 || N > x.size()) {
		N = x.size();
	}
	for(int k=0; k<N; k++) {
		double y;
		for(int n=0; n<x.size(); n++) {
			y += x[n] * cos((M_PI/N) * (n+0.5) * k);
		}
		output.push_back(sqrt(2)/x.size()*y);
	}
	return output;
}

vector<bool> nt_to_binary(string seq) {
	vector<bool> bin;
	for(char c : seq) {
		switch(c) {
			case 'A':
				bin.push_back(0); bin.push_back(0); break;
			case 'C':
				bin.push_back(0); bin.push_back(1); break;
			case 'G':
				bin.push_back(1); bin.push_back(0); break;
			case 'T':
				bin.push_back(1); bin.push_back(1); break;
			default: cerr << "nucleotide unaccepted in sequence " << c << endl;
		}
	}
	return bin;
}

/// average nucleotide frequency
array<float, 4> extract(string s) {
	int last_pos[] = {0,0,0,0};
	int count[] = {0,0,0,0};
	array<float, 4> freq;
	for(int i=0; i < s.length(); i++) {
		switch(s.at(i)) {
			case 'A': count[0]++; break;
			case 'C': count[1]++; break;
			case 'G': count[2]++; break;
			case 'T': count[3]++; break;
			default: cerr << "nucleotide unaccepted in sequence " << s << endl;
		}
	}
	for(int i=0; i<4; i++) {
		freq[i] = 1.0 * count[i] / s.length();
	}
	return freq;
}

/// decompose a nucleotide sequence into binary channels
array<vector<bool>, 4> decompose(string sequence) {
	array<vector<bool>, 4> channels;
	for(auto it=sequence.begin(); it!=sequence.end(); ++it) {
		auto signal = encode_nt(*it);
		for(int i=0; i<4; i++) {
			channels[i].push_back(signal[i]);
		}
	}
	return channels;
}

int main(const int argc, const char* argv[]) {
	auto channels = decompose(string(argv[1]));
	auto frequencies = extract(string(argv[1]));
	
	for(int i=0; i<4; i++) {
		cerr << "FREQ=" << frequencies[i] << " BINARY_SIGNAL= c(";
		for(auto s : channels[i]) {
			cerr << s << ",";
		}
		cerr << ")" << endl << "dct = c(";
		for(auto s : dct(channels[i])) {
			cerr << s << ",";
		}
		cerr << ")" << endl;
	}

	auto binary = nt_to_binary(string(argv[1]));
	cerr << "binary = c(";
	for(auto s : binary) {
		cerr << s << ",";
	}
	cerr << ")" << endl;
	return 0;
}
