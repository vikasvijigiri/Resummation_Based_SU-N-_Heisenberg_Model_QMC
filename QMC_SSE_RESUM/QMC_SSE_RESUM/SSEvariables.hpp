#ifndef _SSEVARIABLES_HPP_DEFINED_
#define _SSEVARIABLES_HPP_DEFINED_

#include <random>
#include <ctime>
#pragma once

extern int lx, ly, N, Ns, Nb, n1, Lc;
extern double J, Q;
extern int nbins, isteps, iter;
extern double Beta, prob_in, prob_rm;
extern int** bsites;


class SSEvariables {
	public:
      	int** lattice_sites( );
      	void declare_variables(int );
      	void set_temperatures(int );

};

class ran {
  private:
    std::mt19937 mt;
    std::uniform_real_distribution < double > dist;

  public:
    ran(double lower = 0.0, double upper = 1.0): mt(std::time(nullptr)), dist(lower, upper) {}
    double operator()() {
    return dist(mt);
  }
};
#endif
