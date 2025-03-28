// compiled with  g++ CN_search.cpp -lgmp -O3

#include <gmp.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>

int main()
{
  // the below has the following quantities hard-coded:
  // The upper bound of the computation:      lines 33,34
  // The quantity P:                          line 38
  // The quantity L = lambda(P)               lines 44,45
  // The array holding the primes dividing L  lines 46
  // e = v_2( LCM( P-1, L0 ) )                line 54

  // version 1: 401 s, mpz_class version
  // version 2: 376 s, machine words where possible, import to mpz_t for exponentiation
  // version 3: 338 s, all arithmetic in mpz_t
  // so we go with all arithmetic in mpz_t

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  auto t1 = high_resolution_clock::now();

  // set search bound, set as 10^24
  mpz_t bound;
  mpz_init_set_ui( bound, 10 );
  mpz_pow_ui( bound, bound, 24 );

  // set preproduct
  mpz_t P;
  mpz_init_set_ui( P, 515410417841 );

  //mpz_init_set_ui( P, 67491742686737 );

  // set L = lambda(P)
  mpz_t L;
  mpz_init_set_ui( L, 115920 );
  uint64_t L64 = 115920;
  #define NUM_BASES 5
  uint64_t fermat_bases[NUM_BASES] = { 2, 3, 5, 7, 23};

  //mpz_init_set_ui( L, 36756720 );
  //uint64_t L64 = 36756720;
  //uint64_t fermat_bases[7] = { 2, 3, 5, 7, 11, 13, 17 };

  // need power of 2 dividing LCM( P-1, L )
  // for the stronger fermat exponent
  int32_t exp_on_2 = 4;
  int32_t pow_of_2 = ( 1 << exp_on_2 );

  // compute r^* = p^{-1} mod L
  // this is the start of  R = (r^* + kL) w/ k = 0
  mpz_t r_star;
  mpz_init(r_star);
  mpz_invert(r_star, P, L);

  // having computed r_star, we now use uint64_t for this quantity
  uint64_t r_star64;
  mpz_export( &r_star64, 0, 1, sizeof(uint64_t), 0, 0, r_star);

  // This is the start of n = Pr^* + kPL w/ k = 0
  mpz_t n;
  mpz_init(n);
  mpz_mul( n, P, r_star);

  mpz_t( strong_exp );
  mpz_init( strong_exp );

  // common difference for n
  mpz_t PL;
  mpz_init( PL );
  mpz_mul( PL, P, L );

  // will need more bases later
  // use bases from the prime divisors of L
  mpz_t base;
  mpz_init( base );

  mpz_t( gcd_result );
  mpz_init( gcd_result );
  // storage for the result of the exponentiation
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

  while( mpz_cmp( n , bound ) < 0 )
  {
    R_composite_factors.push( r_star64 );
    R_prime_factors.clear();

    int i = 0; //counter for fermat_bases

    do
    {
      // set up strong base:  truncated divsion by 2^e means the exponent holds (n-1)/(2^e)
      mpz_tdiv_q_2exp( strong_exp, n, exp_on_2 );
      mpz_set_ui( base, fermat_bases[i] );
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
          mpz_import (r_factor, 1, 1, sizeof(uint64_t), 0, 0, &temp );

          // check gcd before prime testing
          // result1 holds the algebraic factor assoicated with b^((n-1)/(2^e)) + 1
          mpz_add_ui( result1, result1, 1);
          mpz_gcd( gcd_result, result1, r_factor);

          // check that gcd_result has a nontrivial divisor of r_factor
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
          std::cout << "n = " ;
          gmp_printf( "%Zd", n);
          std::cout << " and R = " << r_star64 << " has " << R_composite_factors.size() << " composite factors and " << R_prime_factors.size() << " prime factors." << std::endl;
        std::cout << "and is a base-" << fermat_bases[i] << " Fermat psp." << std::endl;
      }

      // get next Fermat base
      i++;
      // do it again if
      // the number is a Fermat psp and
      // R_composite queue is not empty
    }
    while( is_fermat_psp && !R_composite_factors.empty() && i < NUM_BASES );

    // empty queue
    while( !R_composite_factors.empty() ){ R_composite_factors.pop(); }

    // move to next candidate in arithmetic progression for n and R
    mpz_add( n, n, PL);
    r_star64 += L64;
  }

  mpz_clear( P );
  mpz_clear( L );
  mpz_clear( r_star );
  mpz_clear( n );
  mpz_clear( strong_exp );
  mpz_clear( PL );
  mpz_clear( base );
  mpz_clear( gcd_result );
  mpz_clear( r_factor );
  mpz_clear( result1 );
  mpz_clear( result2 );
  mpz_clear( bound );

  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_double.count() << "ms\n";



}
