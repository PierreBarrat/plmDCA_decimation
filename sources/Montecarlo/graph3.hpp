#ifndef GRAPH3_HPP
#define GRAPH3_HPP

#include <iostream>
#include <string>
#include "mvector.hpp"
#include "graphs.hpp"
#include "graph1.hpp"
#include "graph3.hpp"

class Graph3
{
	public:
	Graph3(size_t n, size_t q) :
		n(n), q(q), L(xstd::mshape<4>(n, n, q, q)), l(xstd::mshape<2>(n, q)) {}

	Graph3(Graph1 const & g1);

	size_t n, q;
	xstd::mvector<4, double> L;
	xstd::mvector<2, double> l;

	//void randomize(double beta = 1.);

	//void randomize_gauss(double beta = 1.);

	void read(std::istream & is);

	std::ostream & print_distribution(std::ostream & os);

	std::ostream & sample_from_distribution(std::ostream & os, size_t m);

	std::ostream & sample_from_distribution_montecarlo(std::ostream & os, size_t m, size_t mc_iters0, size_t mc_iters, std::string const & out_energies_name);

	std::ostream & print_parameters(std::ostream & os);

	void print_parameters(FILE * of);
};

#endif // GRAPH3_HPP
