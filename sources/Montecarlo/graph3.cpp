#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include "mvector.hpp"
#include "gasdev.hpp"
#include "graph1.hpp"
#include "graph2.hpp"
#include "graph3.hpp"

using namespace std;
using namespace xstd;

extern ostream & log_out;

Graph3::Graph3(Graph1 const & g1) : n(g1.n), q(g1.q), L(mshape<4>(n, n, q, q)), l(mshape<2>(n, q))
{
	log_out << "converting Ls" << endl;
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			/*cerr << "K:" << endl;
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					cerr << g1.K[i][j][yi][yj] << " ";
				}
				cerr << endl;
			}*/
			for (size_t yi = 0; yi < q - 1; ++yi) {
				for (size_t yj = 0; yj < q - 1; ++yj) {
					L[i][j][yi][yj] = g1.K[i][j][yi][yj]
						- g1.K[i][j][yi][q - 1]
						- g1.K[i][j][q - 1][yj]
						+ g1.K[i][j][q - 1][q - 1];
					L[j][i][yj][yi] = L[i][j][yi][yj];
				}
			}
		}
	}
	log_out << "converting Hs" << endl;
#pragma omp parallel for
	for (size_t i = 0; i < n; ++i) {
		vector<double> sumKi(q - 1);
		for (size_t j = 0; j < n; ++j) if (j != i) {
			for (size_t yi = 0; yi < q - 1; ++yi) {
				sumKi[yi] += g1.K[i][j][yi][q - 1] - g1.K[i][j][q - 1][q - 1];
			}
		}
		for (size_t yi = 0; yi < q - 1; ++yi) {
			l[i][yi] = sumKi[yi];
		}
	}
	log_out << "conversion done" << endl;
}

void Graph3::read(istream & is)
{
	log_out << "reading Graph3" << endl;
	string tok;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					is >> tok;
					assert(tok == "J");
					is >> tok;
					assert(atoi(tok.c_str()) == int(i));
					is >> tok;
					assert(atoi(tok.c_str()) == int(j));
					is >> tok;
					assert(atoi(tok.c_str()) == int(yi));
					is >> tok;
					assert(atoi(tok.c_str()) == int(yj));
					is >> tok;
					L[i][j][yi][yj] = atof(tok.c_str());
					L[j][i][yj][yi] = L[i][j][yi][yj];
				}
			}
			for (size_t y = 0; y < q; ++y) {
				assert(fabs(L[i][j][y][q - 1]) < 1e-4);
				assert(fabs(L[i][j][q - 1][y]) < 1e-4);
			}
		}
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t yi = 0; yi < q; ++yi) {
			is >> tok;
			assert(tok == "h");
			is >> tok;
			assert(atoi(tok.c_str()) == int(i));
			is >> tok;
			assert(atoi(tok.c_str()) == int(yi));
			is >> tok;
			l[i][yi] = atof(tok.c_str());
		}
		assert(fabs(l[i][q - 1]) < 1e-4);
	}
	log_out << "done" << endl;
}

ostream & Graph3::print_distribution(ostream & os)
{
	vector<size_t> conf(n);

	double norm = 0;
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			x += l[i][conf[i]];
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += L[i][j][conf[i]][conf[j]];
			}
		}
		norm += exp(x);
		size_t j = 0;
		while (j < n && ++conf[j] == q) {
			conf[j] = 0;
			j++;
		}
		if (j == n) {
			break;
		}
	}
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			x += l[i][conf[i]];
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += L[i][j][conf[i]][conf[j]];
			}
		}
		os << "G3 " << exp(x) / norm << endl;
		size_t j = 0;
		while (j < n && ++conf[j] == q) {
			conf[j] = 0;
			j++;
		}
		if (j == n) {
			break;
		}
	}
	return os;
}


ostream & Graph3::sample_from_distribution(ostream & os, size_t m)
{
	vector<size_t> conf(n);

	size_t num = size_t(round(pow(double(q), double(n))));

	vector<double> cumulative(num + 1);

	size_t c = 1;
	while (true) {
		//for (size_t i = 0; i < n; ++i) { cerr << conf[i] << " "; } cerr << endl;
		double x = 0;
		for (size_t i = 0; i < n; ++i) {
			x += l[i][conf[i]];
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; ++j) {
				x += L[i][j][conf[i]][conf[j]];
			}
		}
		double nnp = exp(x);
		cumulative[c] = cumulative[c - 1];
		cumulative[c++] += nnp;
		//cout << c - 2 << ": ( "; for (size_t i = 0; i < n; ++i) cout << conf[i] << " "; cout << ") : " << exp(x) << " " << cumulative[c - 1] << endl;

		size_t j = 0;
		while (j < n && ++conf[j] == q) {
			conf[j] = 0;
			j++;
		}
		if (j == n) {
			break;
		}
	}
	assert(c == num + 1);

	double norm = cumulative[num];

	size_t s = 0;
	while (s < m) {
		double x = norm * drand48();
		//double x = norm * (double(s) / m);

		size_t i0 = 0;
		size_t i1 = num;
		while (i1 != i0 + 1) {
			size_t i = i0 + (i1 - i0) / 2;
			if (x < cumulative[i]) {
				i1 = i;
			} else {
				i0 = i;
			}
		}
		//cout << "x=" << x << " i0=" << i0 << " : ";
		size_t qq = i0;
		for (size_t i = 0; i < n; ++i) {
			os << qq % q;
			if (i < n - 1) {
				os << " ";
			} else {
				os << endl;
			}
			qq /= q;
		}
		++s;
	}
	return os;
}

ostream & Graph3::sample_from_distribution_montecarlo(ostream & os, size_t m, size_t mc_iters0, size_t mc_iters, string const & out_energies_name)
{
	size_t ts = 0;
	vector<size_t> conf(n);
	for (size_t i = 0; i < n; ++i) {
		conf[i] = size_t(q * drand48());
		assert(conf[i] < q);
		//cerr << conf[i] << " ";
	}
	/*ifstream in_conf("conf.tmp");
	for (size_t i = 0; i < n; ++i) {
		in_conf >> conf[i];
		assert(conf[i] < q);
	}
	in_conf.close();*/

	ofstream out_energies(out_energies_name.c_str());
	log_out << "computing initial energy... ";
	double en = 0.;
	for (size_t i = 0; i < n; ++i) {
		en -= l[i][conf[i]];
		for (size_t j = i + 1; j < n; ++j) {
			en -= L[i][j][conf[i]][conf[j]];
		}
	}
	log_out << "done." << endl;
	out_energies << "-1 " << en << endl;

	log_out << "initialize montecarlo sampling... ";
	double tot_de = 0;
	for (size_t k = 0; k < mc_iters0; ++k) {
		size_t i = size_t(n * drand48());
		size_t dq = 1 + size_t((q - 1) * drand48());

		size_t q0 = conf[i];
		size_t q1 = (q0 + dq) % q;

		double e0 = -l[i][q0];
		for (size_t j = 0; j < n; ++j) if (j != i) {
			e0 -= L[i][j][q0][conf[j]];
		}
		double e1 = -l[i][q1];
		for (size_t j = 0; j < n; ++j) if (j != i) {
			e1 -= L[i][j][q1][conf[j]];
		}
		double de = e1 - e0;
		if ((de < 0) || (drand48() < exp(-de))) {
			//cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
			conf[i] = q1;
			tot_de += de;
		}
	}
	log_out << " [tot_de=" << tot_de << "] done." << endl;
	en += tot_de;
	out_energies << "0 " << en << endl;
	tot_de = 0.;
	for (size_t s = 0; s < m; ++s) {
		for (size_t k = 0; k < mc_iters; ++k) {
			size_t i = size_t(n * drand48());
			size_t dq = 1 + size_t((q - 1) * drand48());

			size_t q0 = conf[i];
			size_t q1 = (q0 + dq) % q;

			double e0 = -l[i][q0];
			for (size_t j = 0; j < n; ++j) if (j != i) {
				e0 -= L[i][j][q0][conf[j]];
			}
			double e1 = -l[i][q1];
			for (size_t j = 0; j < n; ++j) if (j != i) {
				e1 -= L[i][j][q1][conf[j]];
			}
			double de = e1 - e0;
			if ((de < 0) || (drand48() < exp(-de))) {
				//cerr << i << " -> " << q1 << " (" << dq << ")" << endl;
				conf[i] = q1;
				tot_de += de;
			}
		}
		log_out << "\rs=" << ++ts << "/" << m << " de=" << tot_de << "                 ";
		out_energies << s << " " << en + tot_de << endl;
		for (size_t i = 0; i < n; ++i) {
			os << conf[i];
			if (i < n - 1) {
				os << " ";
			} else {
				os << endl;
			}
		}
	}
	log_out << endl;
	out_energies.close();
	return os;
}

ostream & Graph3::print_parameters(ostream & os)
{
	log_out << "printing parameters J" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					os << "J " << i << " " << j << " " << yi << " " << yj << " " << L[i][j][yi][yj] << endl;
				}
			}
		}
	}
	log_out << "printing parameters H" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t yi = 0; yi < q; ++yi) {
			os << "h " << i << " " << yi << " " << l[i][yi] << endl;
		}
	}
	log_out << "done" << endl;
	return os;
}

void Graph3::print_parameters(FILE * of)
{
	log_out << "printing parameters J" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			for (size_t yi = 0; yi < q; ++yi) {
				for (size_t yj = 0; yj < q; ++yj) {
					fprintf(of, "J %lu %lu %lu %lu %g\n", i, j, yi, yj, L[i][j][yi][yj]);
				}
			}
		}
	}
	log_out << "printing parameters H" << endl;
	for (size_t i = 0; i < n; ++i) {
		for (size_t yi = 0; yi < q; ++yi) {
			fprintf(of, "h %lu %lu %g\n", i, yi, l[i][yi]);
		}
	}
	log_out << "done" << endl;
}
