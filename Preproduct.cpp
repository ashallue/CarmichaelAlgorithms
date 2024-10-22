#include "Preproduct.h"
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>

// could consider a constructor that passes in the factors of P and lambda(P)
// does not check that init_LofP is actually lambda(P)
// right now this uses trial division
// intended use case is relatively small initializing preproducts
Preproduct::Preproduct( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound )
{
  	uint64_t temp;
  	temp = init_preproduct;

  	mpz_init( P ) ;
  	mpz_import( P, 1, 1, sizeof(uint64_t), 0, 0, &temp );
  	L = init_LofP;
	append_bound = init_append_bound;

	P_len = 0;

  	// factory by trial division
  	if( init_preproduct != 1 )
  	{
		temp = 1;
	  	while( temp*temp < init_preproduct )
	  	{
		  	temp+= 2;
		  	if( ( init_preproduct % temp ) == 0 )
		  	{
			  	init_preproduct = init_preproduct / temp;
			  	P_primes[ P_len ] = temp;
			  	P_len++;
		  	}
	  	}
	  	// possible remaining factor
	  	if( init_preproduct > 1 )
	  	{
		  	P_primes[ P_len ] = init_preproduct;
		  	P_len++;
	  	}
	}
	
	if( init_LofP == 1 )
  	{ 
		L_len = 0;
   	}
	else
  	{
     //factor L here for L_exponents, L_distinct_primes, and L_len
   	}
  	len_appended_primes = 0;     
}


// admissibility check with no gcd
// if the while loop is not taken, this will execute with less than 10 instructions
// assumes len_appended_primes > 0
bool Preproduct::is_admissible( uint64_t prime_to_append )
{
  // if this works how we want it to, this while loop will not be entered often
  while( prime_to_append > next_inadmissible[0] )
  {
    // add 2*p or 4*p and avoid divisibility be 3.  Flip the state of mod 3 status.
    next_inadmissible[0] += (appended_primes[0] << mod_three_status[0]);
    mod_three_status[0] = (mod_three_status[0] == 1 ) ? 2 : 1;
    int i= 1;
    // put next_inadmissible in sorted order having increased the first element
    while( i < len_appended_primes && next_inadmissible[ i-1 ] < next_inadmissible[ i ] )
    {
      std::swap( next_inadmissible[ i-1 ], next_inadmissible[ i ] );
      std::swap( appended_primes[ i-1 ], appended_primes[ i ] );
      std::swap( mod_three_status[ i-1 ], mod_three_status[ i ] );
      i++;
    }
  }
  // now prime to append is less than or equal to next_inadmissible
  // return true if it is less than; return false if it is equal
  return ( prime_to_append < next_inadmissible[0] ) ;
}


int main(void) {

   Preproduct P0( 143, 60, 13 ); 

   std::cout << "Initializing Preproduct : " << 143 << std::endl;
   

   return 0;
}
