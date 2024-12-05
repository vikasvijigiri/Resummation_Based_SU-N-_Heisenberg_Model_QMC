#include "SSEvariables.hpp"
#include <fstream>

int lx, ly, N, Ns, Nb, n1, Lc;
int nbins, isteps, iter;
int **bsites;
double J, Q;
double Beta, prob_in, prob_rm;

void
SSEvariables::declare_variables(int myrank)
{
	std::ifstream vars("input_params.in");
	vars >> lx;
	vars >> J;
	vars >> Q;
	vars >> nbins;
	vars >> iter;
	vars >> isteps;
	ly = 1;				
	Ns = lx * ly;		 				
	Nb = Ns;
	N  = 2;				
}

void
SSEvariables::set_temperatures(int myrank)
{
	double T = 0.02 + (myrank)*0.02;
	Beta = 1 / T;
	prob_in = (J / N) * (Beta * Nb);
	prob_rm = 1.0/prob_in;
}

int **SSEvariables::lattice_sites()
	{
		// Chain lattice
		bsites = new int * [Nb];
	 	for (int i = 0; i < Nb; i++) bsites[i] = new int[2];
		for (int y1=0; y1<ly; ++y1){
		    for (int x1=0; x1<lx; ++x1){
			int s1 = x1 + y1 * lx;
			int x2 = (x1 + 1) % lx;
			int y2 = y1;
			bsites[s1][0] = std::min(s1,x2 + y2 * lx);				// H bonds,  horizontal above
			bsites[s1][1] = std::max(s1,x2 + y2 * lx);	
		// Square lattice		    
		    /*
		        int s1 = x1+y1*lx;
		        int x2 = (x1+1) % lx;
		        int y2 = y1;
		        bsites[s1][0] = s1;
		        bsites[s1][1] = x2+y2*lx;
		        x2 = x1;
		        y2 = (y1+1) % ly;
		        bsites[s1+Ns][0] = s1;
		        bsites[s1+Ns][1] = x2+y2*lx;
		    */    
		    }
		}        
		return bsites;
	}
