#include "Preproduct.h"
#include <algorithm>
#include <iostream>
#include <queue>
#include <vector>
#include <bit>
#include <bitset>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>
#include <cstddef> 

static_assert(sizeof(unsigned long) == 8, "unsigned long must be 8 bytes");

Preproduct::Preproduct()
{
    mpz_init( P ) ;
    mpz_init( L ) ;
}

Preproduct::~Preproduct()
{
    mpz_clear( P );
    mpz_clear( L );
}

// assumes valid inputs: 
// 1) does not check that init_preproduct is cylic 
// 2) does not check that init_LofP is actually CarmichaelLambda( init_preproduct )
// intended use is initializing from precomputation which only generates valid inputs
// uses trial division because the intended use case is relatively small initializing preproducts
// could consider another version that passes in arrays of factors
// precomputation.cpp would need to be re-written
void Preproduct::initializing( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound )
{
    mpz_set_ui( P, init_preproduct );
    mpz_set_ui( L, init_LofP );
    append_bound = init_append_bound;
    P_len = 0;

       
    uint64_t temp;
    // set primes array for P
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
    
    // set primes and exponent arrays for L
    if( init_LofP == 1 ) { L_len = 0; }
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
                L_distinct_primes[ L_len ] = temp;
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
void Preproduct::appending( Preproduct PP, primes_stuff p )
{
    mpz_mul_ui( P, PP.P, p.prime );
    P_len = PP.P_len + 1;
    std::copy( PP.P_primes,PP.P_primes + PP.P_len, P_primes );
    P_primes[ P_len ] =  p.prime;   
    append_bound = p.prime;
    
    // initialize L to be PP.L and increase only when new factors are seen
    // L = PP.L;
    // if( (PP.L) % (p.prime - 1) == 0 )

    mpz_set( L, PP.L );
    if( mpz_divisible_ui_p( L, ( p.prime - 1 ) )  )
    {
        std::copy( PP.L_distinct_primes,PP.L_distinct_primes + PP.L_len, L_distinct_primes );
        std::copy( PP.L_exponents, PP.L_exponents + PP.L_len, L_exponents );
        L_len = PP.L_len;
    }
    else
    {
        //merge L and p-1
        int i = 0; //counter for PP
        int j = 0; //counter for p
        L_len = 0; //counter for P
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
                        // need to import PP.L_distinct_primes[i] before multplying ?
                        mpz_mul_ui( L, L, PP.L_distinct_primes[i] );
                        //L *= PP.L_distinct_primes[i];
                    }
                }
                i++; j++; L_len++;
            }
            else if( PP.L_distinct_primes[i] < p.pm1_distinct_primes[j] )
            {
                L_distinct_primes[ L_len ] = PP.L_distinct_primes[i];
                L_exponents[ L_len ] = PP.L_exponents[i];
                i++; L_len++;
            }
            else
            {
                L_distinct_primes[ L_len ] = p.pm1_distinct_primes[j];
                L_exponents[ L_len ] = p.pm1_exponents[j];
                for( int L_update = 0; L_update < p.pm1_exponents[j]; L_update++ )
                {
                    mpz_mul_ui( L, L, p.pm1_distinct_primes[j] );
                    //L *= p.pm1_distinct_primes[j];
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
                mpz_mul_ui( L, L, p.pm1_distinct_primes[j] );
                // L *= p.pm1_distinct_primes[j];
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
        temp_mod_3 = 1; 
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

// some analysis could be done to minimize mpz_init calls
// other optimizations in the if( is_fermat_psp ) branch
// data structure choice? queue right now
// check modular exponentation prior to computing gcd
// note that these optimizations effect a minority of the computation
void Preproduct::CN_search( uint64_t bound_on_R )
{
    // there are two arithmetic progressions associated with n = P*R
    // letting r^* = P^{-1} mod L where 0 < r^* < L
    // 1) R = r^* + kL - common difference of L
    // 2) n = PR = Pr^* + kPL - common difference of PL

    
    mpz_t r_star;
    mpz_init( r_star );
    // using r_star as a temporary variable for set-up
    
    mpz_set( r_star, P );
    // it now holds the same value as P
    mpz_sub_ui( r_star, r_star, 1);
    // it now holds P-1
    
    // need power of 2 dividing LCM( P-1, L )
    // for the stronger fermat exponent
    // we could dynamically choose for each n
    // but we choose the largest power of 2 that works for all n
    // does this matter?  if so, this needs to be moved to the primary while loop
    int32_t exp_on_2 = std::min( L_exponents[ 0 ], (uint16_t) mpz_scan1( r_star, 0) );
    int32_t pow_of_2 = ( 1 << exp_on_2 );

    // compute r^* = p^{-1} mod L
    // this is the start of  R = (r^* + kL) w/ k = 0
    // r_star is no longer being used as a temporary variable
    // it now it holds the correct value
    mpz_invert(r_star, P, L);

    // having computed r_star, we now use uint64_t for this quantity
    uint64_t r_star64;
    mpz_export( &r_star64, 0, 1, sizeof(uint64_t), 0, 0, r_star);

    uint64_t L64;
    mpz_export( &L64, 0, 1, sizeof(uint64_t), 0, 0, L);
    
    // This is the start of n = Pr^* + kPL w/ k = 0
    // so n = Pr^*
    mpz_t n;
    mpz_init( n );
    mpz_mul( n, P, r_star);

    mpz_t strong_exp;
    mpz_init( strong_exp );

    // common difference for n
    mpz_t PL;
    mpz_init( PL );
    mpz_mul( PL, P, L );

    // will need more bases later
    // use bases from the prime divisors of L
    mpz_t base;
    mpz_init( base );

    // storage for the gcd result
    mpz_t gcd_result;
    mpz_init( gcd_result );
    
    // storage for the result of the exponentiation
    // result1 will hold the stronger test
    // result2 will hold the Fermat test
    mpz_t result1;
    mpz_init( result1 );
    mpz_t result2;
    mpz_init( result2 );

    mpz_t r_factor;
    mpz_init( r_factor );

    bool is_fermat_psp;

    std::queue<uint64_t> R_composite_factors;
    std::vector<uint64_t> R_prime_factors;

    uint64_t temp;

    while( r_star64 <= bound_on_R )
    {
      R_composite_factors.push( r_star64 );
      R_prime_factors.clear();

      int i = 0; //counter for fermat_bases

      do
      {
        // set up strong base:  truncated divsion by 2^e means the exponent holds (n-1)/(2^e)
        mpz_tdiv_q_2exp( strong_exp, n, exp_on_2 );
        // we use prime divisors of L as the Fermat bases
        mpz_set_ui( base, L_distinct_primes[ i ] );
        mpz_powm( result1,  base,  strong_exp, n); // b^( (n-1)/(2^e) )
        mpz_powm_ui( result2,  result1, pow_of_2, n); // b^( (n-1)/(2^e)) )^(2^e) = b^(n-1)

        is_fermat_psp = ( mpz_cmp_si( result2, 1 ) == 0 );

        // this conditional is not expected to be entered
        // so the do-while loop is not expected to be invoked
        // most numbers are not Fermat pseudoprimes
        if( is_fermat_psp )
        {
          int start_size = R_composite_factors.size();
          // use a for loop to go through all factors that are currently in the queue
          for( int j = 0; j < start_size; j++ )
          {
            // get element out of queue and put into mpz_t
            // first time through, this is just r_factor will have the value of r_star
            temp = R_composite_factors.front();
            R_composite_factors.pop();
            
            // should be mpz_import ?
            mpz_set_ui( r_factor, temp );
              
            // check gcd before prime testing
            // result1 holds the algebraic factor assoicated with b^((n-1)/2^e) + 1
            mpz_add_ui( result1, result1, 1);
            mpz_gcd( gcd_result, result1, r_factor);

            // check that gcd_result has a nontrivial divisor of r_factor
            // could probably be a check on result1 = +/- 1 mod n
            // before computing the gcd
            if( mpz_cmp(gcd_result, r_factor) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
            {
              // will need to add a check about a lower bound on these divisors
              mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result);
              ( mpz_probab_prime_p( gcd_result, 0 ) == 0 ) ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
              mpz_divexact(gcd_result, r_factor, gcd_result );
              mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result);
              ( mpz_probab_prime_p( gcd_result, 0 ) == 0 ) ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
            }
            else // r_factor was not factored, so it is prime or composite
            {
              ( mpz_probab_prime_p( r_factor, 0 ) == 0 ) ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
            }
          }
          // if R_composite is empty, check n is CN *here*
         
          // 
          gmp_printf ("n = %Zd", n);
          std::cout << " and R = " << r_star64 << " has " << R_composite_factors.size() << " composite factors and " << R_prime_factors.size() << " prime factors." << std::endl;
          std::cout << "and is a base-" << L_distinct_primes[i] << " Fermat psp." << std::endl;
        }

        // get next Fermat base
        i++;
        // do it again if
        // the number is a Fermat psp and
        // R_composite queue is not empty
        // if i == L_len, we should probably output or factor directly 
            // could be some strange multi-base Fermat pseudoprime - very rare?
        // room for improvement here
      }
      while( is_fermat_psp && !R_composite_factors.empty() && i < L_len );

      // empty queue
      while( !R_composite_factors.empty() ){ R_composite_factors.pop(); }

      // move to next candidate in arithmetic progression for n and R
      mpz_add( n, n, PL);
      r_star64 += L64;
    }

    mpz_clear( r_star );
    mpz_clear( n );
    mpz_clear( strong_exp );
    mpz_clear( PL );
    mpz_clear( base );
    mpz_clear( gcd_result );
    mpz_clear( r_factor );
    mpz_clear( result1 );
    mpz_clear( result2 );
    
}

bool Preproduct::appending_is_CN( std::vector< uint64_t >&  primes_to_append )
{
    mpz_t P_temp;
    mpz_t L_temp;
    mpz_t gcd_result;

    mpz_init_set( P_temp, P );
    mpz_init_set( L_temp, L );
    mpz_init_set_ui( gcd_result, 1);
    
    bool return_val = true;

    // not meaningful optimization availalbe
    // break the loop(s) when return_val is detected to be false
    for( auto app_prime : primes_to_append )
    {
        for( int i = 0; i < P_len; i++ )
        {
            return_val = ( return_val && ( app_prime % P_primes[i] != 1 ) );
        }
        // update the preproduct by this prime
        mpz_mul_ui( P_temp, P_temp, app_prime );
        uint64_t temp = app_prime - 1;
        // compute LCM( L, p-1 ) = (L / gcd( L, p-1 ) )*(p-1)
        mpz_gcd_ui( gcd_result, L_temp, temp );
        mpz_divexact( L_temp, L_temp, gcd_result );
        mpz_mul_ui( L_temp, L_temp, temp );
    }

    mpz_sub_ui( P_temp, P_temp, 1);
    return_val = return_val && mpz_divisible_p( P_temp, L_temp );
        
    mpz_clear( P_temp );
    mpz_clear( L_temp );
    mpz_clear( gcd_result );
    
    return return_val;
    
}

//
bool Preproduct::is_CN( )
{
    bool return_val;
   
    // subtract 1 from P
    mpz_sub_ui( P, P, 1);
    
    // check if L exactly divides P-1
    return_val = mpz_divisible_p( P  , L );
    
    // return P to its correct value
    mpz_add_ui( P, P, 1);
    
    return return_val;
}



int main(void) {
    
    Preproduct P0;
    P0.initializing( 6682828353, 2289560, 13 );
    
    std::cout << "LP for this is " << sizeof( unsigned long int ) << std::endl;
    
    std::cout << "Initializing P : " ;
    gmp_printf ("%Zd = ", P0.P );
    
    for( int i = 0; i < ( P0.P_len - 1 ); i++)
    {
        std::cout << P0.P_primes[i] << " * "  ;       
    }
    std::cout << P0.P_primes[P0.P_len - 1 ] << std::endl ;  
    
    std::cout << "Initializing Lambda : " ;
    gmp_printf ("%Zd = ", P0.L );
    
    for( int i = 0; i < ( P0.L_len - 1 ); i++)
    {
        std::cout << P0.L_distinct_primes[i] << " ^ "  << P0.L_exponents[ i ] << " * "  ;       
    }
    std::cout << P0.L_distinct_primes[P0.L_len - 1] << " ^ "  << P0.L_exponents[ P0.L_len - 1 ] << std::endl ;  
    
    P0.CN_search(149637241475922);
    
    
    return 0;
}
