#ifndef _SSEUPDATES_HPP_DEFINED_
#define _SSEUPDATES_HPP_DEFINED_
#include <tuple>
class SSEupdates {
      	public: 
      	SSEupdates(){}
	int  *str, *frst, *last, *X, *lspn; 


      	void initialize( );
	void mcstep( );
	void checkl( );
	
	void update( );
	void looper( );		
	std::tuple <int, int, int, int, int>  calculate_dnl(int , int , int, bool );	   	
	void update_links(int, int, int, int, int, int, int, bool );
	
};
#endif

