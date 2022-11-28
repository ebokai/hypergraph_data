// METROPOLIS DATA GENERATOR FOR 
// LARGE SYSTEMS WITH GENERAL INTERACTIONS


// specify list of interactions (in integer format)
// specify list of parameters 

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <bitset>
#include <cmath>

using namespace std;

const int n = 20;

map<uint32_t, float> read_interactions(string &int_file);
void metropolis_sampling(int &N, int &multiplier, map<uint32_t, float> &interactions, string &out_file);


int main(int argc, char **argv){

	int N;
	int multiplier = 5;
	string interaction_file;
	string output_file;

	sscanf(argv[1], "%d", &N);
	interaction_file = argv[2];
	output_file = argv[3];

	cout << "generating: " << N << " data-points" << endl;
	cout << "loading interactions from: " << interaction_file << endl;
	cout << "writing data to: " << output_file << ".dat" << endl;
	
	map<uint32_t, float> interactions = read_interactions(interaction_file);

	map<uint32_t, float>::iterator it;

	// for (it = interactions.begin(); it != interactions.end(); it++){
	// 	cout << it->first << " " << it->second << endl;
	// }

	metropolis_sampling(N, multiplier, interactions, output_file);

	return 0;

}

map<uint32_t, float> read_interactions(string &int_file){

	// read interactions and parameters from file

	map<uint32_t, float> interactions;

	ifstream myfile(int_file);
	string line, subline;

	while(getline(myfile, line)){
		istringstream ss(line);
		string token;
		float linedata[2];
		int i = 0;

		while(getline(ss, token, ';')){
			linedata[i] = stof(token);
			i++;

		interactions[(uint32_t) linedata[0]] = linedata[1];
			
		}
	}

	myfile.close();

	return interactions;

} 

void metropolis_sampling(int &N, int &multiplier, map<uint32_t, float> &interactions, string &out_file){

	srand(time(NULL));

	uint32_t state = 0;
	uint32_t new_state, op;
	// vector<uint32_t> data;
	map<uint32_t, float>::iterator it;
	float e0, e1, g, delta, u;

	int eval0, eval1;
	int steps = N * n * multiplier;

	ofstream data;
	string fname = out_file + ".dat";
	data.open(fname);


	for (int t = 0; t < steps; t++){

		// flip spin - XOR with 2^i
		uint32_t i = rand() / (RAND_MAX/n);

		uint32_t spin = pow(2,i);
		new_state = (state ^ spin);

		e0 = 0;
		e1 = 1;

		for (it = interactions.begin(); it != interactions.end(); it++){

			op = it->first;
			g = it->second;

			eval0 = 1 - 2 * (bitset<n>(op & state).count() % 2);
			eval1 = 1 - 2 * (bitset<n>(op & new_state).count() % 2);

			e0 += g * eval0;
			e1 += g * eval1;
		}

		delta = e0 - e1; // this depends on convention 

		u = static_cast <float> (rand()) / static_cast <float> (RAND_MAX); 

		if (delta <= 0){
			state = new_state;
		}
		else if (exp(-delta) > u){
			state = new_state;
		}

		if (t % (n * multiplier) == 0){
			data << bitset<n>(state);
			data << "\n";
		}	

	}

	data.close();
}