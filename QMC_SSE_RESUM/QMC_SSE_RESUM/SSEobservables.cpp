#include "SSEvariables.hpp"
#include "SSEupdates.hpp"
#include "SSEobservables.hpp"
#include <iostream>

void
SSEobservables::observables( )
{
	//double m = 0;
	//int i,j,b,s;
	/*for (int i=0; i < Ns; ++i){
            m += lattice[i][mb[s][0]].S()*pow(-1,((i-1) % lx)+(i-1)/lx);
            }
        m/=2;
	*/
	double mag = 0.;
	double mag_abs = 0.;
	double mag_square = 0.;
	int jj[2] = { 0 };
	/*
        for (int p=0; p < Lc; ++p){
            if (upd.str[0][p] == 0) continue;
            else if (upd.str[0][p] == 2){
                 b=upd.str[1][p]-1;
                 s=upd.str[2][p]-1;
                 i=bsites[b][0];
                 j=bsites[b][1];
                	//std::cout << i << "  " << j << "    " << b << std::endl; 
                 lattice[i][mb[s][0]].flip();
                 lattice[j][mb[s][1]].flip();
                 jj[int((b-1)/Ns)]+=lattice[j].S();
                 m+=2*lattice[i][mb[s][0]].S()*pow(-1,((i-1) % lx)+(i-1)/lx);
                 }
            mag+=float(m);
            mag_abs+=float(abs(m));
            mag_square+=pow(float(m),2);
        }    
        if (n1 == 0){
           mag=(pow(mag,2)+mag_square)/(float(n1)*float(n1+1));
           mag_abs/=n1;
           mag_square/=n1;
	}
	else{
           mag_abs=float(abs(m));
           mag_square=pow(float(m),2);
           mag=mag_square;
        } 
        */
	enrg1 += float(n1);
	enrg2 += pow(float(n1), 2);
	amag_abs += mag_abs;
	amag_square += mag_square;
	asusc = asusc + mag;
	stiff += 0.5 *(pow(float(jj[0]), 2) + pow(float(jj[1]), 2));
	//ususc+=pow(float(sum(lattice[:].S())/2),2);
}

void
SSEobservables::binning_data(int bins, int iters)
{

	enrg1 /= iters;
	enrg2 /= iters;
	amag_abs /= iters;
	amag_square /= iters;
	asusc /= iters;
	stiff /= iters;
	ususc /= iters;

	enrg2 = (enrg2 - enrg1 *(enrg1 + 1.)) / Ns;
	enrg1 = -enrg1 / (Beta *Ns);
	amag_abs = amag_abs / Ns;
	amag_square = amag_square / Ns;
	asusc = Beta *asusc / Ns;
	ususc = Beta *ususc / Ns;
	stiff = stiff / (Beta *Ns);

	//std::cout << "energ " << enrg1 << std::endl; 
	data1[0] += enrg1;
	data1[1] += enrg2;
	data1[2] += amag_abs;
	data1[3] += amag_square;
	data1[4] += asusc;
	data1[5] += stiff;
	data1[6] += ususc;

	data2[0] += pow(enrg1, 2);
	data2[1] += pow(enrg2, 2);
	data2[2] += pow(amag_abs, 2);
	data2[3] += pow(amag_square, 2);
	data2[4] += pow(asusc, 2);
	data2[5] += pow(stiff, 2);
	data2[6] += pow(ususc, 2);

	for (int i = 0; i < 7; ++i)
	{
		wdata1[i] = data1[i] / bins;
		wdata2[i] = data2[i] / bins;
		wdata2[i] = sqrt(abs(wdata2[i] - pow(wdata1[i], 2)) / bins);
	}

	enrg1 = 0.;
	enrg2 = 0.;
	amag_abs = 0.;
	amag_square = 0.;
	asusc = 0.;
	stiff = 0.;
	ususc = 0.;
}

void
SSEobservables::Initiate_observables()
{
	data1 = new double[6];
	data2 = new double[6];
	wdata1 = new double[6];
	wdata2 = new double[6];
	for (int i = 0; i < 7; ++i)
	{
		data1[i] = 0.;
		data2[i] = 0.;
	}
}
