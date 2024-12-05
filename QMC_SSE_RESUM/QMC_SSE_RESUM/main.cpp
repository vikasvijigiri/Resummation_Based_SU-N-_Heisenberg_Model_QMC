#include "SSEvariables.hpp"
#include "SSEupdates.hpp"			//Headers
#include "SSEobservables.hpp"
#include "writeresults.hpp"
#include "mpi.h"
#include <iostream>

static ran rann;
using namespace std;
int main()
{
	SSEvariables vars;	// Creating all necessary instances
	SSEupdates update;
	SSEobservables obs;
	writeresults write;

						
	MPI::Status status;	// MPI part
	MPI::Init();
	int myrank = MPI::COMM_WORLD.Get_rank();
	int numprocs = MPI::COMM_WORLD.Get_size();	

	vars.declare_variables(myrank);
	vars.set_temperatures(myrank);

	double t1=0, t2=0;
	if (myrank == 0)
	{
		t1 = MPI::Wtime();
	}


	bsites = vars.lattice_sites();
	update.initialize();

  	update.looper();
	for (int i = 1; i <= isteps; ++i)
	{
		// Equilibration 
		update.mcstep();
		update.checkl();
	} 
	obs.Initiate_observables();
	
	for (int j = 1; j <= nbins; ++j)
	{
		// Measurements
		for (int i = 1; i <= iter; ++i)
		{
			update.mcstep();
			obs.observables();
		}
		obs.binning_data(nbins, iter);
	}
	write.output(myrank, 1 / Beta, obs.wdata1[0], obs.wdata1[1], obs.wdata1[2], obs.wdata1[3], obs.wdata1[4], obs.wdata1[5], obs.wdata1[6], numprocs);
	if (myrank == 0)
	{
		t2 = MPI::Wtime();
		std::cout << "Total time taken by the prcoesses is, " << (t2 - t1) / 60 << ", minutes" << std::endl;
	}
	MPI::Finalize();
	return 0;
}
