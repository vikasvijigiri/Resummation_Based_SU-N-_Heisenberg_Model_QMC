#include "SSEvariables.hpp"
#include "writeresults.hpp"
#include <iostream>
#include <fstream>
#include "mpi.h"
using namespace std;

//***************************************************************************************
void writeresults::output(int myrank, double h_local, double enrg_local, double mag_local, double m_local, double p_local, double J1_local, double op_local, double Ic_local, int numprocs)
{
	double *h_global    = NULL;
	double *mag_global  = NULL;
	double *enrg_global = NULL;
	double *m_global    = NULL;
	double *p_global    = NULL;
	double *J1_global   = NULL; 
	double *Ic_global   = NULL; 
	double *op_global   = NULL; 
	
	h_global    = (double*) malloc(numprocs* sizeof(double));
	mag_global  = (double*) malloc(numprocs* sizeof(double));
	enrg_global = (double*) malloc(numprocs* sizeof(double));
	m_global    = (double*) malloc(numprocs* sizeof(double));
	p_global    = (double*) malloc(numprocs* sizeof(double));
	J1_global   = (double*) malloc(numprocs* sizeof(double));
	Ic_global   = (double*) malloc(numprocs* sizeof(double));
	op_global   = (double*) malloc(numprocs* sizeof(double));

	MPI_Gather(&h_local,    1, MPI::DOUBLE, h_global,    1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&mag_local,  1, MPI::DOUBLE, mag_global,  1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&enrg_local, 1, MPI::DOUBLE, enrg_global, 1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&m_local,    1, MPI::DOUBLE, m_global,    1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&p_local,    1, MPI::DOUBLE, p_global,    1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&J1_local,   1, MPI::DOUBLE, J1_global,   1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&Ic_local,   1, MPI::DOUBLE, Ic_global,   1, MPI::DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&op_local,   1, MPI::DOUBLE, op_global,   1, MPI::DOUBLE, 0, MPI::COMM_WORLD);

	if (myrank == 0)
	{
		std::cout << "Energy per spin at T = " << 1 / Beta << "  is, " << enrg_global[0] << "  " << enrg_local << std::endl;
		std::string resultfile;
		resultfile += "../QMC_SSE_RESUM/results/data_" + to_string_with_precision(J, 2) + "_" + to_string_with_precision(Q, 2) + "_" +
		 to_string_with_precision(lx, 2) + ".dat";

		std::ofstream file;
		file.open(resultfile);

		//if (!file);
		for (int i = 0; i < numprocs; ++i)
		{
			file << setprecision(6) << fixed;	//set some precision for nice format
			file << setw(15) << h_global[i];
			file << setw(15) << mag_global[i];
			file << setprecision(5) << fixed;	     
			file << setw(15) << enrg_global[i];
			file << setw(15) << m_global[i];
			file << setw(15) << p_global[i];
			file << setw(15) << J1_global[i];
			file << setw(15) << op_global[i];
			file << setw(15) << Ic_global[i];
			file << '\n';
		}

		file.close();
	}
	delete h_global;
	delete mag_global;
	delete enrg_global;
	delete m_global;
	delete p_global;
	delete J1_global;
	delete Ic_global;
	delete op_global;
}

//***************************************************************************************
template < typename A>
	std::string writeresults::to_string_with_precision(A a_value, int n)
	{
		std::ostringstream out;
		out << std::setprecision(n) << a_value;
		return out.str();
	}
//***************************************************************************************
