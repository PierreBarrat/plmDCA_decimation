#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <unistd.h>
#include "graphs.hpp"
#include "graph1.hpp"
#include "graph2.hpp"
#include "graph3.hpp"

using namespace std;

void usage(char * command, ostream & os)
{

	os << "Usage: " << command << " [options]" << endl;
	os << endl;
	os << "Montecarlo." << endl;

	os << "  options:" << endl;
	os << "    -n <n>    number of nodes" << endl;
	os << "    -q <q>    number of states" << endl;
	os << "    -m <m>    number of samples" << endl;
	os << "    -s <s>    random seed" << endl;
	os << "    -t <t>    montecarlo iterations (default=100000)" << endl;
	os << "    -T <t>    montecarlo initial iterations (default=1000000)" << endl;
	os << "    -0        use gauge0" << endl;
	os << "    -u <str>  extra suffix for output files" << endl;
	os << "    -?        this help" << endl;
}

int main(int argc, char ** argv)
{
	int c;
	size_t n = 4;
	size_t q = 3;

	size_t m = 10;

	long int seed = 1;

	size_t mc_iters0 = 1000000;
	size_t mc_iters = 100000;

	bool gauge0 = false;

	string suffix("");

	while ((c = getopt(argc, argv, "n:q:m:s:b:t:T:0?u:")) != -1) {
		switch (c) {
			case 'n':
				n = atoi(optarg);
				break;
			case 'q':
				q = atoi(optarg);
				break;
			case 'm':
				m = atoi(optarg);
				break;
			case 's':
				seed = atol(optarg);
				break;
			case 'T':
				mc_iters0 = atoi(optarg);
				break;
			case 't':
				mc_iters = atoi(optarg);
				break;
			case '0':
				gauge0 = true;
				break;
			case 'u':
				suffix = string(optarg);
				break;
			case '?':
				usage(argv[0], cout);
				exit(EXIT_SUCCESS);
				break;
			default:
				cerr << "Use -? for help." << endl;
				exit(EXIT_FAILURE);
		}
	}
	srand48(seed);

	string out_samples_name;

	stringstream ss;
	ss << "out_samples_montecarlo_n" << n << "_q" << q << "_s" << seed << "_m" << m << "_T" << mc_iters0 << "_t" << mc_iters << suffix << ".txt";
	ss >> out_samples_name;

	ofstream out_samples(out_samples_name.c_str());

	string out_energies_name;

	stringstream se;
	se << "out_energies_montecarlo_n" << n << "_q" << q << "_s" << seed << "_m" << m << "_T" << mc_iters0 << "_t" << mc_iters << suffix << ".txt";
	se >> out_energies_name;

	out_samples << m << " " << n << " " << q << endl;

	if (not gauge0) {
		Graph2 g2(n, q);
		g2.read(cin);

		//cout << "-------------" << endl;
		//g2.print_distribution(cout);

		g2.sample_from_distribution_montecarlo(out_samples, m, mc_iters0, mc_iters, out_energies_name);
	} else {
		Graph3 g3(n, q);
		g3.read(cin);

		//cout << "-------------" << endl;
		//g3.print_distribution(cout);

		g3.sample_from_distribution_montecarlo(out_samples, m, mc_iters0, mc_iters, out_energies_name);
	}
	out_samples.close();

	return EXIT_SUCCESS;
}
