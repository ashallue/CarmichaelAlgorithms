#include "Preproduct.h"
#include <algorithm>
#include <iostream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>


// assumes valid inputs: 
// 1) does not check that P is cylic 
// 2) does not check that init_LofP is actually CarmichaelLambda(P)
// uses trial division because the intended use case is relatively small initializing preproducts
// could consider a constructor that passes in the factors of P and lambda(P)
Preproduct::Preproduct( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound )
{
  	uint64_t temp;
  	temp = init_preproduct;

  	mpz_init( P ) ;
  	mpz_import( P, 1, 1, sizeof(uint64_t), 0, 0, &temp );
  	L = init_LofP;
	append_bound = init_append_bound;

	P_len = 0;

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
		L_len = 1;
	    L_distinct_primes[0] = 2;
	    L_exponents[0] = __builtin_ctzl(init_LofP );
		init_LofP = ( init_LofP >>  __builtin_ctzl( init_LofP ) );
		// the prime 2 is now taken care of.  Account for odd primes:
		temp = 3;
		while( init_LofP != 1 )
		{	
			if( init_LofP % temp == 0 )
			{
				L_exponents[ L_len ] = 0;
				L_distinct_primes[L_len] = temp;
				while( init_LofP % temp == 0 )
				{
					L_exponents[ L_len ]++;
					init_LofP /= temp;
				}
				L_len++;
			}	
			temp += 2;
		}
	}
  	len_appended_primes = 0;     
}

// assumes prime_stuff is valid and admissible to PP
Preproduct::Preproduct( Preproduct PP, primes_stuff p )
{
	mpz_init( P ) ;
	mpz_mul_ui( P, PP.P, p.prime );
	// preproduct now has the correct value

	P_len = PP.P_len + 1;
	std::copy( PP.P_primes,PP.P_primes + PP.P_len, P_primes );
	P_primes[ P_len ] =  p.prime;
	// primes array now has the correct primes
	
	append_bound = p.prime;
	// append bound is updated
	
	// initialize L to be PP.L and increase only when new factors are seen
	L = PP.L;
	
	if( (PP.L) % (p.prime - 1) == 0 )
	{
		std::copy( PP.L_distinct_primes,PP.L_distinct_primes + PP.L_len, L_distinct_primes );
		std::copy( PP.L_exponents, PP.L_exponents + PP.L_len, L_exponents );
		L_len = P.L_len;
	}
	else
	{
		//merge L and p-1
		int i = 0; //counter for PP
		int j = 0; //counter for p
		L_len = 0; // counter for P
		while( i < PP.L_len && j < p.pm1_len )
		{
			// check if they have the same prime
			if( PP.L_distinct_primes[i] == p.pm1_distinct_primes[j] )
			{
				L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
				L_exponents[ L_len ] = std::max( PP.L_exponents[i], p.pm1_exponents[j] );
				//we need to update L if  L_exp < pm1_exp
				if( PP.L_exponents[i] <  p.pm1_exponents[j])
				{
					for( int L_update = PP.L_exponents[i]; L_update < p.pm1_exponents[j]; L_update++ )
					{
						L *= PP.L_distinct_primes[i];
					}
				}
				i++; j++; L_len++;
			}
			else if( PP.L_distinct_primes[i] < p.pm1_distinct_primes[j] )
			{
				L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
				L_exponents[ L_len ] = PP.L_exponents[i];
				// no update to L required
				i++; L_len++;
			}
			else
			{
				L_distinct_primes[ L_len ] = p.pm1_distinct_primes[j];
				L_exponents[ L_len ] = p.pm1_exponents[j];
				for( int L_update = 0; L_update < p.pm1_exponents[j]; L_update++ )
				{
					L *= p.pm1_distinct_primes[j];
				}
				j++; L_len++;
			}
		}
		while( i < PP.L_len )
		{
			L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
			L_exponents[ L_len ] = PP.L_exponents[i];
			i++; L_len++;
		}	
		while( j < p.pm1_len  )
		{
			L_distinct_primes[ L_len ] = p.pm1_distinct_primes[j];
			L_exponents[ L_len ] = p.pm1_exponents[j];
			for( int L_update = 0; L_update < p.pm1_exponents[j]; L_update++ )
			{
				L *= p.pm1_distinct_primes[j];
			}
			j++; L_len++;
		}
	}
	
	//set appended prime info for further admissibility checks
	len_appended_primes = PP.len_appended_primes + 1;
	// the next inadmissible prime to p.prime is 2*p.prime + 1 or 4*p.prime + 1
	// case is chosen to avoid divisibility by 3
	uint64_t temp_next_inad = 2*p.prime + 1;
	uint16_t temp_mod_3 = 2;
	if( p.prime % 3 == 1 )
	{
		temp_next_inad += 2*p.prime;
		temp_mod_3 = 1; //next time add 2^1 * p
	}
	
	if( len_appended_primes == 1 )
	{
	    next_inadmissible[0] = temp_next_inad; 
	    mod_three_status[0] = temp_mod_3;  
	    appended_primes[0] = p.prime; 
	}
	else
	{
		int i = 0;
		while( PP.next_inadmissible[i] < temp_next_inad )
		{
		    next_inadmissible[i] = PP.next_inadmissible[i];
		    mod_three_status[i] = PP.mod_three_status[i];
		    appended_primes[i] = PP.appended_primes[i];
			i++;
		}
	    next_inadmissible[i] = temp_next_inad; 
	    mod_three_status[i] = temp_mod_3;  
	    appended_primes[i] = p.prime; 
		while( i < PP.len_appended_primes )
		{
		    next_inadmissible[i+1] = PP.next_inadmissible[i];
		    mod_three_status[i+1] = PP.mod_three_status[i];
		    appended_primes[i+1] = PP.appended_primes[i];
			i++;
		}

	}
	

}

Preproduct::~Preproduct()
{
	mpz_clear( P ); 
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

   std::cout << "Initializing Lambda : " << P0.L << " =  " ;
   
   for( int i = 0; i < ( P0.L_len - 1 ) ; i++ )
   {
	   std::cout << P0.L_distinct_primes[i] << " ^ "  << P0.L_exponents[ i ] << " * "  ;		
   }
   std::cout << P0.L_distinct_primes[P0.L_len - 1] << " ^ "  << P0.L_exponents[ P0.L_len - 1 ] << std::endl ;	

   return 0;
}
