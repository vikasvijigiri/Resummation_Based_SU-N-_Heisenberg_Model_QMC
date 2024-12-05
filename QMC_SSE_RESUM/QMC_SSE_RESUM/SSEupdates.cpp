#include "SSEvariables.hpp"
#include "SSEupdates.hpp"
#include "fastmath.hpp"
#include <iostream>
#include <tuple>

static ran rann;
using namespace std;

//********************************************************************************************
//********************************************************************************************

void
SSEupdates::mcstep ()
{
  update ();
}

void
SSEupdates::checkl()
{
	if (int(n1 + n1 / 3) > Lc)
	{
		int Lc_new = int(n1 + n1 / 3);
		int *str_tmp = new int[Lc_new];
		int *X_tmp = new int[4 *Lc_new];
		
		for (int i = 0; i < Lc_new; ++i)
		{
			if (i < Lc)
			{
				str_tmp[i] = str[i];
				X_tmp[4 *i] = X[4 *i];
				X_tmp[4 *i + 1] = X[4 *i + 1];
				X_tmp[4 *i + 2] = X[4 *i + 2];
				X_tmp[4 *i + 3] = X[4 *i + 3];
				/*
				if (l1 < -1 or l2 < -1 or l3 < -1 or l4 < -1 ){
				cout << l1 << "  " << l2 << "  " << l3 << "  " << l4 << endl;
				printf(" X error  ");
				exit(0);  
				}
				else if (l1 > 4*Lc or l2 > 4*Lc  or l3 > 4*Lc  or l4 > 4*Lc  ){
				cout << l1 << "  " << l2 << "  " << l3 << "  " << l4 << endl;
				printf(" X error  ");
				exit(0);  
				}
				*/
			}
			else
			{
				str_tmp[i] = 0;
				X_tmp[4 *i] = -1;
				X_tmp[4 *i + 1] = -1;
				X_tmp[4 *i + 2] = -1;
				X_tmp[4 *i + 3] = -1;
			}
		}
		//cout << "couldn't come till here c " << endl;
		delete[] str;
		delete[] X;

		//cout << "couldn't come till here d " << endl;

		Lc = Lc_new;
		str = new int[Lc];
		X = new int[4 *Lc];

		for (int i = 0; i < Lc; ++i)
		{
			str[i] = str_tmp[i];
			X[4 *i] = X_tmp[4 *i];
			X[4 *i + 1] = X_tmp[4 *i + 1];
			X[4 *i + 2] = X_tmp[4 *i + 2];
			X[4 *i + 3] = X_tmp[4 *i + 3];
		}
		delete[] str_tmp;
		delete[] X_tmp;
		//cout << "New adjusted cutoff length is " << Lc << endl;
	}
}



//********************************************************************************************

void
SSEupdates::initialize ()
{
  Lc = 20;
  n1 = 5;
  int b;
  
  
  frst = new int [Ns];
  last = new int [Ns];
  X    = new int [4*Lc];  
  str  = new int [Lc];

  for (int j = 0; j < Lc; j++) str[j] = 0;
  
  b=1;
  str[0]=b+1;

  b=3;
  str[2]=b+1;
  
  b=2;
  str[10]=b+1;

  b=1;
  str[14]=b+1;
  
  b=2;
  str[17]=b+1;
}



std::tuple <int, int, int, int, int> SSEupdates::calculate_dnl (int s1, int s2, int i, bool insert) {

  int dnl=0;
  
  int vi, v[2], v1[2], w[2], pos_type[2], l0, l2, id=0, iu=0;

  v[0] = frst[s1];
  v[1] = frst[s2];
  
  w[0] = last[s1];
  w[1] = last[s2];
  
  v1[0] = v[0];
  v1[1] = v[1];
  
  l0 = 4 * i;
  l2 = l0 + 2;
 
  /*
  if(v[0] != X[X[v[0]]] and v[0] != -1) {cout << "start links v " << v[0] << "  " << X[v[0]] << "  " << X[X[v[0]]] << endl; exit(0);  }
  if(v[1] != X[X[v[1]]] and v[1] != -1) {cout << "start links v " << v[1] << "  " << X[v[1]] << "  " << X[X[v[1]]] << endl; exit(0);  }
  if(w[0] != X[X[w[0]]] and w[0] != -1) {cout << "start links w " << w[0] << "  " << X[w[0]] << "  " << X[X[w[0]]] << endl; exit(0);  } 
  if(w[1] != X[X[w[1]]] and w[1] != -1) {cout << "start links w " << w[1] << "  " << X[w[1]] << "  " << X[X[w[1]]] << endl; exit(0);  }       
  
  if (v[0] == v[1] and v[0] != -1) { printf(" error in the frst legs, they are same! "); exit(0); }
  if (w[0] == w[1] and w[0] != -1) { printf(" error in the last legs, they are same! "); exit(0); }  
  */
  
  // ********************************** Getting the leg vertex of the site. **************************************
  if (insert) {
     for (int j = 0; j < 2; j++)
  	{
  		if (v[j] != -1)
  		{
  		        id = int(v[j]/4);
        		iu = int(w[j]/4);  
  			if ( i < id and i < iu )	// link ABOVE to the current one
  			{	
  				v1[j] = v[j];
  				pos_type[j] = 0;
  			}
  			else if ( i > id and i < iu )	// link SANDWITCH to the current one.
  			{
  				v1[j] = v[j];
  				//cout << "upp legs are   " << v1[j] << endl;  
  				while (1)
  				{
  					v1[j] += 2;
  					v1[j] %= 4*Lc; 
  					//cout << "upp legs are   " << v1[j] << endl;  					
  					v1[j] = X[v1[j]];
  					//cout << "upp legs are   " << v1[j] << endl;  					
  					//if (v1[j] == v[j]) break;
  					if (int(v1[j]/4) > i)
  					{
  						v1[j] = X[v1[j]];
  						//v1[j] -= 2;
  						break;
  					}
  				}
  				/*
  				if ( (abs(4 * i - v1[j]) > abs(4 * i - v1[j] - 2)) and (i < int(Lc / 2)) )
  				{
  					v1[j] += 2;  					
  					v1[j] %= 4*Lc;
  				}
  				*/
  				pos_type[j] = 1;
  			}
  			else if ( i > id and i > iu )	// link BELOW to the current one
  			{
  				v1[j] = w[j];
  				pos_type[j] = 2;
  			}
  			else 
  			{
    				//cout << "legs are   " << j << "  " << i << "  " << id << "  " << iu << endl;   		 			
  	        		printf(" Error in getting vertex leg! ");
				exit(0);  			
  			}
  		}
  		else
  		{
  				pos_type[j] = 3;	
  		}
        }
  }
  else
  {
     for (int j = 0; j < 2; j++)
  	{
         	id = int(X[(l0+j)]/4);
        	iu = int(X[(l2+j)]/4);  	
		if ( i == id and  i == iu ){
			pos_type[j] = 3; // None 
		}  		
		else if ( (i > id and i < iu) or (i < id and i > iu) ){
			pos_type[j] = 1; // Sand
		}
		else if ( i > id and i > iu ){
			pos_type[j] = 2; // Below
		}	
		else if ( i < id and i < iu ){
			pos_type[j] = 0; // Above		
		}
  		else {
   		//cout << "legs are   " << i << "  " << v1[0] << "  " << v1[1] << "  " << id << "  " << iu << endl;   		
  	        	printf(" Error in computing vertex leg! ");
			exit(0);  			
  		}				
		v1[j] = l0 + 2*j;
	}  		
  }
  
  //cout << "legs are   " << v1[0] << "  " << v1[1] << "  " << pos_type[0] << "  " << pos_type[1] << "  " << id << "  " << iu << endl;  
  /*
  if ((v1[0] < 0  and pos_type[0] < 3) or (v1[1] < 0 and pos_type[1] < 3)){
	        printf(" Error in the pos_type index! ");
		exit(0);  
  }
  */

  // ******************************************* Check SAME loop or NOT. ***********************************************
  vi = v1[0];
  //int count = 0;
  if (v1[0] != -1 and v1[1] != -1)
  {
    		//cout << "dwn legs are   " << vi << endl;  	
  	while (1)
  	{
  		//count += 1;
  		vi ^= 1;
  		//cout << "dwn legs are   " << vi << endl;  	  		
  		if (vi == v1[1])	// same loop 
  		{
  			dnl = 1;
  			break;
  		}		
  		vi = X[vi];
  		//cout << "dwn legs are   " << vi << endl;  	  		
  		if (vi == v1[1])	// same loop
  		{
  			dnl = 1;
  			break;
  		}	 
  		if (vi == v1[0])	// diff loop
  		{
  			dnl = -1;
  			break;
  		}
  		/*
  		if (count > 200){
	        printf(" Stuck here ");
		exit(0);  		
  		}
  		*/	
  	}
  }
  else	// Link NONE above or below to the current one.
  {
  	for (int j = 0; j < 2; j++)
  	{
  		if (v1[j] < 0) 		// diff loop
  		{
  			v1[j] = 4 * i + j;
  		}	
  	}
  	dnl = -1;	
  }
  return {dnl, v1[0], v1[1], pos_type[0], pos_type[1]};		
}


void
SSEupdates::update_links (int s1, int s2, int i, int vl1, int vl2, int pos_type1, int pos_type2, bool insert)
{
  int s[2], v1[2], vu[2], vd[2], pos_type[2];
  int l0, l2;
  
  s[0] = s1;
  s[1] = s2;
  
  v1[0] = vl1;
  v1[1] = vl2;

  pos_type[0]=pos_type1;
  pos_type[1]=pos_type2; 
    
  l0=4*i;
  l2=4*i+2;
  
  for (int j = 0; j < 2; j++){
	  vd[j] = v1[j];
	  if (pos_type[j] == 0)	    // link above.  
	  {
	    if (insert){
	        vu[j] = X[vd[j]];	    

	        X[l0+j]=vu[j];
	        X[vu[j]]=l0+j;
	        
	        X[l2+j]=vd[j];
	        X[vd[j]]=l2+j;
	    
	        frst[s[j]]=l0+j;
	    //cout << "links are above ins  "  << vd[j] << "  " << vu[j] << "  " << X[vd[j]] << "  " << X[vu[j]] << "  " << X[X[vd[j]]] << "  " << X[X[vu[j]]] 
	    // << endl;     	    	        	        
	    }
	    else {
		vd[0] = X[l0 + j];
		vd[1] = X[l2 + j];

		frst[s[j]]=vd[1];

		X[vd[0]] = vd[1];
		X[vd[1]] = vd[0];
				  
		X[l0 + j] = -1;
		X[l2 + j] = -1;
	    //cout << "links are above rem  "  << vd[0] << "  " << vd[1] << "  " << X[vd[0]] << "  " << X[vd[1]] << "  " << X[X[vd[0]]] << "  " << X[X[vd[1]]] 
	    //<< endl;     	    	        	        
	    }
	    //cout << "links are above    " << X[29] << "  " << X[31] << "  " << frst[s1] << endl;     	    	        	    	    
	  }
	  else if (pos_type[j] == 1) 	    // link Sandwitching.
	  {
	    if (insert){	  
	        vu[j] = X[vd[j]];
	    
	        X[l0+j]=vd[j];
	        X[vd[j]]=l0+j;
	    
	        X[l2+j]=vu[j];
	        X[vu[j]]=l2+j;
	    //cout << "links are sandw ins  "  << vd[j] << "  " << vu[j] << "  " << X[vd[j]] << "  " << X[vu[j]] << "  " << X[X[vd[j]]] << "  " << X[X[vu[j]]] 
	    // << endl;     	    	        
	    }
	    else {
		vd[0] = X[l0 + j];
		vd[1] = X[l2 + j];

		X[vd[0]] = vd[1];
		X[vd[1]] = vd[0];
		  
		X[l0 + j] = -1;
		X[l2 + j] = -1;
	    //cout << "links are sandw rem  "  << vd[0] << "  " << vd[1] << "  " << X[vd[0]] << "  " << X[vd[1]] << "  " << X[X[vd[0]]] << "  " << X[X[vd[1]]] 
	    //<< endl;      	    	        	        
	    } 
	  }  
	  else if (pos_type[j] == 2)	    // link below.
	  {
	    if (insert){
	    	vu[j] = X[vd[j]];
	    
	    	X[l0+j]=vd[j];
	    	X[vd[j]]=l0+j;

	    	last[s[j]]=l2+j;
	    
	    	X[l2+j]=vu[j];
	    	X[vu[j]]=l2+j;
	    //cout << "links are below ins  "  << vd[j] << "  " << vu[j] << "  " << X[vd[j]] << "  " << X[vu[j]] << "  " << X[X[vd[j]]] << "  " << X[X[vu[j]]] 
	    // << endl;      	    	        
	    }
	    else {
		vd[0] = X[l0 + j];
		vd[1] = X[l2 + j];

	    	last[s[j]] = vd[0];

		X[vd[0]] = vd[1];
		X[vd[1]] = vd[0];

		X[l0 + j] = -1;
		X[l2 + j] = -1;	
	    //cout << "links are below rem  "  << vd[0] << "  " << vd[1] << "  " << X[vd[0]] << "  " << X[vd[1]] << "  " << X[X[vd[0]]] << "  " << X[X[vd[1]]] 
	    //<< endl;      	    	        	        
	    }
	    //cout << "links are below     " << X[16] << "   " << X[18] << "  " << last[s2] << endl;     	    	    	        	  	  
	  }
	  else if (pos_type[j] == 3) // None above or below.  
	  {
	    if (insert){
	    	X[l0 + j]  = l2 + j;
	    	X[l2 + j]  = l0 + j;
	    	frst[s[j]] = l0 + j;
	    	last[s[j]] = l2 + j;
	    }
	    else {
	        vl1 = X[l0 + j];
	        vl2 = X[l2 + j];
		X[l0 + j]  = -1;
		X[l2 + j]  = -1;		    
	    	last[s[j]] = -1;
	    	frst[s[j]] = -1;
	    }
	   //cout << "links are None     " << X[0] << "   " << X[1] << "   " << X[2] << "   " << X[3] << endl;     	    	    	        	 
	  } 
	  else {
	        printf(" Error in the pos_type index! ");
		exit(0);	  
	  }	  
  }  
}



void
SSEupdates::update() {
  int b, o, s1, s2;
  double p;

  for (int i = 0; i < Lc; ++i) {
    o = str[i];
    if (o == 0){
      b = int(rann() * Nb);     
      s1 = bsites[b][0];
      s2 = bsites[b][1];
      auto [dnl, v1, v2, pos_type1, pos_type2] = calculate_dnl(s1, s2, i, true);
      p = prob_in * pow(N, dnl);
      if (p > float(Lc-n1) or p > rann()*(Lc-n1)){
        str[i] = b + 1;
        n1 += 1; 
        update_links(s1, s2, i, v1, v2, pos_type1, pos_type2, true);
      }
    } else {
      b  = str[i] - 1;
      s1 = bsites[b][0];
      s2 = bsites[b][1];
      auto [dnl, v1, v2, pos_type1, pos_type2] = calculate_dnl(s1, s2, i, false);
      p = prob_rm * pow(N, dnl);
      if (p*(Lc-n1+1) > 1.0  or p*(Lc-n1+1) > rann()){
        str[i] = 0;
        n1 -=  1; 
        update_links(s1, s2, i, v1, v2, pos_type1, pos_type2, false);
      }
    }
  }
}

void
SSEupdates::looper() {
  /*
  ---------------- Linked list-storage-----------
  */
  std::fill_n(frst, Ns, -1);
  std::fill_n(last, Ns, -1);
  std::fill_n(X, 4 * Lc, -1);

  int o, v0, v1, v2, s1, s2, b;
  for (int i = 0; i < Lc; ++i) {
    v0 = 4*i;
    o  = str[i];
    if (o > 0) {
      b = str[i] - 1;
      s1 = bsites[b][0];
      s2 = bsites[b][1];

      v1 = last[s1];
      v2 = last[s2];
      if (v1 != -1) {
        X[v1] = v0;
        X[v0] = v1;
      } else {
        frst[s1] = v0;
      }
      if (v2 != -1) {
        X[v2] = v0 + 1;
        X[v0 + 1] = v2;
      } else {
        frst[s2] = v0 + 1;
      }
      last[s1] = v0 + 2;
      last[s2] = v0 + 3;
    }
  }
  /*
  ---------------- PBC loops --------------------------!
  */
  for (int k = 0; k < Ns; ++k) {
    v1 = frst[k];
    if (v1 != -1) {
      v2 = last[k];
      X[v2] = v1;
      X[v1] = v2;
    }
  }
  }
