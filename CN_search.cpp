// compiled with  g++ CN_search.cpp -lgmpxx -lgmp -O3

#include <gmpxx.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <queue>
#include <vector>

int main()
{

  // version 1: 401 s, mpz_class version
  // version 2: 376 s, machine words where possible, import to mpz_t for exponentiation
  // version 3: 338 s, all arithmetic in mpz_t
  // so we go with all arithmetic in mpz_t

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  auto t1 = high_resolution_clock::now();

  // set search bound, set as 10^23
  mpz_t bound;
  mpz_init_set_ui( bound, 10 );
  mpz_pow_ui( bound, bound, 24 );

  // set preproduct
  mpz_t P;
  mpz_init_set_ui( P, 6682828353 );

  // set L = lambda(P)
  mpz_t L;
  mpz_init_set_ui( L, 2289560 );

  // primes dividing L will serve as Fermat bases
  uint64_t fermat_bases[6] = { 2, 5, 7, 13, 17, 37 };


  // compute r^* = p^{-1} mod L
  // this is the start of  R = (r^* + kL) w/ k = 0
  mpz_t r_star;
  mpz_init(r_star);
  mpz_invert(r_star, P, L);

  // This is the start of n = Pr^* + kPL w/ k = 0
  mpz_t n;
  mpz_init(n);
  mpz_mul( n, P, r_star);

  // exponent is (n - 1)/2
  // accomplished by using trucated division of n/2
  mpz_t( EJ_exp );
  mpz_init( EJ_exp );
  mpz_tdiv_q_2exp( EJ_exp, n, 1 );

  // common difference for n
  mpz_t PL;
  mpz_init( PL );
  mpz_mul( PL, P, L );

  // will need more bases later
  // use bases from the prime divisors of L
  mpz_t base;
  mpz_init_set_ui( base, fermat_bases[0] );

  // use
  mpz_t other_base;
  mpz_init( other_base );

  mpz_t( gcd_result );
  mpz_init( gcd_result );
  // storage for the result of the exponentiation
  mpz_t result1;
  mpz_init( result1 );
  mpz_t result2;
  mpz_init( result2 );

  bool is_fermat_psp;

  std::queue<uint64_t> R_composite_factors;
  std::vector<uint64_t> R_prime_factors;

  uint64_t temp;

  while( mpz_cmp( n , bound ) < 0 )
  {
    // using base 2:  result1 has EJ test and result2 has Fermat test
    mpz_powm( result1,  base,  EJ_exp, n);
    mpz_powm_ui( result2,  result1, 2, n);

    is_fermat_psp = ( mpz_cmp_si( result2, 1 ) == 0 );

    if( is_fermat_psp )
    {
      // check gcd before prime testing
      // result1 holds the algebraic factor assoicated with b^((n-1)/2) + 1
      mpz_add_ui (result1, result1, 1);
      mpz_gcd (gcd_result, result1, r_star);

      // check that gcd_result has a nontrivial divisor of r_star
      if( mpz_cmp(gcd_result, r_star) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
      {
        // will need to add a check about a lower bound on these divisors
        mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result);
        ( mpz_probab_prime_p( gcd_result, 0 ) == 0 ) ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
        mpz_divexact(gcd_result, r_star, gcd_result );
        mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result);
        ( mpz_probab_prime_p( gcd_result, 0 ) == 0 ) ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
      }
      else // r_star was not factored, so it is prime or composite
      {
        mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, r_star);
        ( mpz_probab_prime_p( r_star, 0 ) == 0 ) ? R_composite_factors.push( temp ) : R_prime_factors.push_back( temp );
      }
      // CN check here if R_composite_factors is empty.
    }

    // we want to enter this loop if Fermat test has passed and R_composite queue is not empty
    while( is_fermat_psp && !R_composite_factors.empty() )
    {

      // For each R_i in R_compsite attempt to factor with gcd with EJ result1
      // use gcd( R_i, result1 + 1 ) or gcd( R_i, result1 - 1)
      // nontrivial creates two factors:  check both for primality
      // put result either into prime or rcf_temp.
      // need to be mindful of the minimum bound that the precomputation requires prime factors to satisfy
      // could have an "abort" here upin finding a factor that is too small or... see ahead

      // if R_composite is empty
      //    check CN with Korselt
      //    we could incorporate the requirement on prime factors in R /w/r/t k here, too.
      // else
      //    get next valid base and try again
      //    mpz_powm( result1,  other_base,  EJ_exp, n);
      //    mpz_powm_ui( result2,  result1,  2, n);

     std::cout << "n = " << n << " and R = " << r_star << " and R is composite." <<  std::endl;
    }

    // move to next candidate in arithmetic progression for n and R
    // update exponent for next round
    mpz_add( n, n, PL);
    mpz_add( r_star, r_star, L);
    mpz_tdiv_q_2exp( EJ_exp, n, 1 );

  }

  mpz_clear( P );
  mpz_clear( L );
  mpz_clear( r_star );
  mpz_clear( n );
  mpz_clear( EJ_exp );
  mpz_clear( PL );
  mpz_clear( base );
  mpz_clear( other_base );
  mpz_clear( gcd_result );
  mpz_clear( result1 );
  mpz_clear( result2 );
  mpz_clear( bound );

  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_double.count() << "ms\n";



}
