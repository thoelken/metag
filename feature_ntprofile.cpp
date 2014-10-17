#include <string>
#include <iostream>
#include <array>
#include <vector>

using namespace std;

array<int, 4> encode_nt(char nt) {
	array<int, 4> encoding = {0,0,0,0};
	switch(nt) {
		case 'A': encoding[0] = 1; break;
		case 'C': encoding[1] = 1; break;
		case 'G': encoding[2] = 1; break;
		case 'T': encoding[3] = 1; break;
		default: cerr << "nucleotide unaccepted in sequence " << nt << endl;
	}
	return encoding;
}

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

array<vector<int>, 4> decompose(string sequence) {
	array<vector<int>, 4> channels;
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
		cerr << "FREQ=" << frequencies[i] << " BINARY_SIGNAL=";
		for(auto s : channels[i]) {
			cerr << s;
		}
		cerr << endl;
	}
	return 0;
}