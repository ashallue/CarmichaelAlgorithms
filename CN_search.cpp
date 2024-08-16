// compiled with  g++ CN_search.cpp -lgmpxx -lgmp -O3

#include <gmpxx.h>
#include <iostream>
#include <chrono>

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
  mpz_pow_ui( bound, bound, 23 );

  // set preproduct
  mpz_t P;
  mpz_init_set_ui( P, 6682828353 );

  // set L = lambda(P)
  mpz_t L;
  mpz_init_set_ui( L, 2289560 );

  // compute r^* = p^{-1} mod L
  mpz_t r_star;
  mpz_init(r_star);
  mpz_invert(r_star, P, L);

  // keep track of two arithmetic progressions
  // R = r^* + kL, L is as above
  // n = P( r^* + k*L) = Pr + kPL
  mpz_t n;
  mpz_init(n);
  mpz_mul( n, P, r_star);

  // exponent is n - 1
  mpz_t( EJ_exp );
  mpz_init( EJ_exp );
  // truncated division by 2 acocmplishes (n-1)/2
  mpz_tdiv_q_2exp( EJ_exp, n, 1 );

  // common difference for n
  mpz_t PL;
  mpz_init( PL );
  mpz_mul( PL, P, L );

  // will need more bases later
  mpz_t base;
  mpz_init_set_ui( base, 2 );

  // storage for the result of the exponentiation
  mpz_t result1;
  mpz_init( result1 );
  mpz_t result2;
  mpz_init( result2 );

  while( mpz_cmp( n , bound ) < 0 )
  {
    // result1 has Euler-Jacobi test
    mpz_powm( result1,  base,  EJ_exp, n);

    // result2 has Fermat test
    mpz_powm_ui( result2,  result1,  2, n);

    // check if we have a Fermat pseudoprime
    if( mpz_cmp_si( result2, 1 ) == 0 )
    {
      // need gcd( R, result1 + 1 ) or gcd( R, result1 - 1)


      std::cout << "n = " << n << " and R = " << r_star << std::endl;
    }

    // next candidate in arithmetic progression
    mpz_add( n, n, PL);

    // truncated division by 2 acocmplishes (n-1)/2
    mpz_tdiv_q_2exp( EJ_exp, n, 1 );

    // and the next R in arithmetic progression
    mpz_add( r_star, r_star, L);
  }

  mpz_clear( P );
  mpz_clear( L );
  mpz_clear( r_star );
  mpz_clear( n );
  mpz_clear( EJ_exp );
  mpz_clear( PL );
  mpz_clear( base );
  mpz_clear( result1 );
  mpz_clear( result2 );
  mpz_clear( bound );

  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_double.count() << "ms\n";



}
