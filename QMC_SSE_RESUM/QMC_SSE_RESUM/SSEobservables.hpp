#ifndef _SSEOBSERVABLES_HPP_DEFINED_
#define _SSEOBSERVABLES_HPP_DEFINED_
#pragma once
class SSEobservables {
  	public:
      	double enrg1, enrg2, amag_abs, amag_square, asusc, stiff, ususc;
      	double *data1, *data2, *wdata1, *wdata2;
      	void observables( );
      	void binning_data(int , int );
      	void Initiate_observables();
};
#endif
