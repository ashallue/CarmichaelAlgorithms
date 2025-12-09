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



  // set search bound, set as 10^24
  mpz_t bound;
  mpz_init_set_ui( bound, 10 );
  mpz_pow_ui( bound, bound, 25 );

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
    std::cout << "r_star = " << r_star64 << std::endl;
    
    
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

    uint32_t f_psp = 0;
    uint32_t mod_exp_count = 0;
    
    mpz_set_ui( base, 2 );
    
    auto t1 = high_resolution_clock::now();
    
    while( mpz_cmp( n , bound ) < 0 )
    {
        mod_exp_count++;
        mpz_powm( result1,  base,  n, n);

        if( mpz_cmp_ui( result1, 2 ) == 0 )
        {
            f_psp++;
        }

        mpz_add( n, n, PL);
    }
    auto t2 = high_resolution_clock::now();
    std::cout << "I found checked " << mod_exp_count << " numbers and found " << f_psp << " base-2 Fermat psps." << std::endl;

    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << " it took " << ms_double.count() << "ms\n";
    
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





}
