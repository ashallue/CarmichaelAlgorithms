#include "Preproduct.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <queue>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>
#include <cstddef>
#include <boost/dynamic_bitset.hpp>

static_assert(sizeof(unsigned long) == 8, "unsigned long must be 8 bytes.  needed for mpz's unsigned longs to take 64 bit inputs in various calls.  LP64 model needed ");


Preproduct::Preproduct()
{
    mpz_init( P );
    mpz_init( L );
    mpz_init_set_ui( BOUND, 10 );
    mpz_pow_ui( BOUND, BOUND, 24 );
    

    
}

Preproduct::~Preproduct()
{
    mpz_clear( P );
    mpz_clear( L );
    mpz_clear( BOUND );
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
    uint64_t temp;
        
    //set info for P:
    mpz_set_ui( P, init_preproduct );
    
    // clear vector and fill
    P_primes.clear();
    if( init_preproduct != 1 )
    {
        temp = 1;
        while( temp*temp < init_preproduct )
        {
            temp+= 2;
            if( ( init_preproduct % temp ) == 0 )
            {
                init_preproduct = init_preproduct / temp;
                P_primes.push_back(temp);
            }
        }
        // possible remaining factor
        if( init_preproduct > 1 )
        {
            P_primes.push_back( init_preproduct );
        }
    }
    
    //set info for L:
    mpz_set_ui( L, init_LofP );
    // clear vector and fill
    L_primes.clear();
    
    // take care of 2 then look for odd primes
    L_primes.push_back( 2 );
    init_LofP = ( init_LofP >>  __builtin_ctzl( init_LofP ) );
    // the prime 2 is now taken care of.  Account for odd primes:
    temp = 3;
    while( init_LofP != 1 )
    {
        if( init_LofP % temp == 0 )
        {
            L_primes.push_back( temp );
            while( init_LofP % temp == 0 )
            {
                init_LofP /= temp;
            }
        }
        temp += 2;
    }
    
    // set append_bound
    append_bound = init_append_bound;

}


// assumes prime_stuff is valid and admissible to PP
void Preproduct::appending( Preproduct PP, primes_stuff p )
{
    // compute new P value:
    mpz_mul_ui( P, PP.P, p.prime );
    // set the vector
    P_primes.clear();
    // copy the old vector over
    P_primes = PP.P_primes;
    // add the new prime
    P_primes.push_back( p.prime );
        
    // compute L
    mpz_lcm_ui( L, PP.L, p.prime - 1 );
    // set L's vector
    L_primes.clear();
    std::set_union( (PP.L_primes).begin(), (PP.L_primes).end(), (p.pm1_distinct_primes).begin(), (p.pm1_distinct_primes).end(), std::back_inserter( L_primes ) );
    
    append_bound = p.prime;
    
}


bool Preproduct::is_admissible_modchecks( uint64_t prime_to_append )
{
    bool return_val = true;
    
    // can be done with std::all_of or std::any_of or std::none_of
    
    uint16_t i = 0;
    while( return_val && i < P_primes.size() )
    {
        return_val = return_val && (1 != ( prime_to_append % P_primes[i] ) );
        i++;
    }
 
    return return_val;
}


void Preproduct::CN_search(  )
{
    const uint32_t cache_bound = 150'000;
        
    mpz_t r_star;
    mpz_init( r_star );
    mpz_invert(r_star, P, L);
    
    mpz_t base2;
    mpz_t base3;
    mpz_init_set_ui( base2, 2 );
    mpz_init_set_ui( base3, 3 );
    
    mpz_t fermat_result;
    mpz_init( fermat_result );
    
    mpz_t n;
    mpz_init(n);

    // n will be created in an arithmetic progression with P*L*(possibly lifted small primes)
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
    
    mpz_t mpz_prime;
    mpz_init( mpz_prime );
    
    mpz_t small_prime;
    mpz_init( small_prime );
    
    mpz_t lifted_L;
    mpz_init_set( lifted_L, L );

    
    // position i has the truth value of the statement "(2i + 3) is prime"
    std::bitset<256> small_primes{"0010010100010100010000010110100010010100110000110000100010100010010100100100010010110000100100000010110100000010000110100110010010010000110010110100000100000110110000110010010100100110000110010100000010110110100010010100110100110010010110100110010110110111"};
    // fix append_bound to be consistent with the bitset
    // we probably need to also bound this by cmp_bound
    // it's not really a sieve if p > cmp_bound
    uint64_t sieve_index_bound = std::min( (append_bound - 3)/2, (uint64_t) 256);
    
  
    // remove primes dividing L from the bitset
    // i initialized to 1 so that the prime 2 is skipped
    uint16_t i = 1;
    while( i < L_primes.size() &&  L_primes[i] < 512 )
    {
        small_primes [ ( L_primes[i] - 3) /2 ] = 0;
        i++;
    }
    
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    mpz_cdiv_q( cmp_bound, BOUND, PL);
   
    mpz_t R;
    mpz_init( R );
    
    // store the primes that we append
    std::vector< uint16_t > primes_lifting_L;
    
    // L_lift will store the product of primes and gives the count of spokes on the wheel
    uint64_t L_lift= 1;
    uint16_t prime_index = 0;
    
    while( mpz_cmp_ui( cmp_bound, cache_bound ) > 0  )
    {
        if( small_primes[ prime_index ] )
        {
            uint16_t p = 2*prime_index + 3;
            L_lift *= p;
            primes_lifting_L.push_back( p );
            mpz_mul_ui( PL, PL, p );
            mpz_mul_ui( lifted_L, lifted_L, p );
            mpz_cdiv_q( cmp_bound, BOUND, PL);
        }
        prime_index++;
    }

    uint32_t cmp_bound32 = mpz_get_ui( cmp_bound );
    boost::dynamic_bitset<> spoke_sieve( cmp_bound32 );
        
    // This is the wheel and a creates r_star + m*L
    // The wheel is r_star + m*L
    for( uint64_t m = 0; m < L_lift; m++ )
    {
        bool enter_loop = true;
        // checks to make sure that r_star + m*L is not divisible by lifted_primes
        for( auto p : primes_lifting_L )
        {
            enter_loop = ( enter_loop && ( mpz_divisible_ui_p( r_star, p ) == 0 ) );
        }
        
        // This is a spoke.  It sieves on k-values numbers of the form (r_star + m*L) + k*L_lift*L
        if( enter_loop )
        {
            //sieve a spoke
            spoke_sieve.reset();
            uint16_t sieve_prime_index = prime_index;
            while( sieve_prime_index < sieve_index_bound )
            {
                // r_star + k*(L_lift*L) = 0 mod p
                // implies k = -r*( L_lift*L)^{-1} mod p
               if( small_primes[ sieve_prime_index ] )
                {
                    uint32_t p = 2*sieve_prime_index + 3;
                    mpz_set_ui( small_prime, p );
                    // n is being used as a temporary variable in the following 5 lines
                    mpz_invert( n, lifted_L, small_prime );     // n has (L_lift*L)^{-1} mod p
                    mpz_neg( n, n);                             // n has -(L)^{-1} mod p
                    mpz_mul( n, n, r_star );                    // muliply by r_star
                    mpz_mod( n, n, small_prime );               // reduce modulo p
                    uint64_t k = mpz_get_ui( n );

                    while( k < cmp_bound32 )
                    {
                        spoke_sieve[k] = 1;
                        k += p;
                    }
                }
                sieve_prime_index++;
            }
            
            // spoke has been sieved, so only do modular exponentiations on valid places
            // n is initialized correctly here
            mpz_mul( n, P, r_star);
            uint32_t k = 0;
            while( mpz_cmp( n , BOUND ) < 0 )
            {
                if( spoke_sieve[k] == 0 )
                {
                    mpz_powm( fermat_result,  base2,  n, n); // 2^n mod n
                    if( mpz_cmp( fermat_result, base2 ) == 0 )  // check if 2 = 2^n mod n
                    {
                        mpz_powm( fermat_result,  base3,  n, n); // 3^n mod n
                        if( mpz_cmp( fermat_result, base3 ) == 0 )  // check if 3 = 3^n mod n
                        {
                            mpz_divexact( R, n, P);
                            CN_factorization( n, R );
                            // gmp_printf( "n = %Zd = %Zd * %Zd is a base-2 and base-3 Fermat psp. \n", n, P, R);
                        }
                    }
                }
                k++;
                mpz_add( n, n, PL);
            }
        }
        mpz_add( r_star, r_star, L);
    }
       
    mpz_clear( R );
    mpz_clear( small_prime );
    mpz_clear( cmp_bound );
    mpz_clear( r_star );
    mpz_clear( n );
    mpz_clear( PL );
    mpz_clear( base2 );
    mpz_clear( base3 );
    mpz_clear( fermat_result );
    mpz_clear( lifted_L );
        
}



bool Preproduct::appending_is_CN( std::vector< uint64_t >&  primes_to_append )
{

    mpz_t P_temp;
    mpz_t L_temp;

    mpz_init_set( P_temp, P );
    mpz_init_set( L_temp, L );

    // "update" P and L for each prime
    for( auto app_prime : primes_to_append )
    {
        // update the preproduct by this prime
        mpz_mul_ui( P_temp, P_temp, app_prime );
        uint64_t temp = app_prime - 1;
        mpz_lcm_ui( L_temp, L_temp, temp );
    }

    mpz_sub_ui( P_temp, P_temp, 1);
    bool return_val =  mpz_divisible_p( P_temp, L_temp );

    mpz_clear( P_temp );
    mpz_clear( L_temp );

    return return_val;
}


// Check that lambda(P) divides (P-1)
// consider changing this to int-type return matching how gmp returns
// and have the same return standard as gmp
bool Preproduct::is_CN( )
{
    bool return_val;
    mpz_sub_ui( P, P, 1);
    return_val = mpz_divisible_p( P, L );
    mpz_add_ui( P, P, 1);
    return return_val;
}


// compare with:
// https://github.com/ashallue/tabulate_car/blob/master/LargePreproduct.cpp#L439C1-L500C2
// see section 5.3 of ANTS 2024 work
// assume that B/PL is "big"
// if the prime to be found is of the form r_star + k*L for some relatively small value of k we should probably invoke a call to CN_search
// the below needs to be fixed so that it is bounded
void Preproduct::completing_with_exactly_one_prime()
{
    // set up scaled problem:
    mpz_t R;
    mpz_init( R );
    
    mpz_t script_R;
    mpz_init( script_R );
    mpz_invert( script_R, P, L );
    mpz_sub_ui( script_R, script_R, 1 );
    
    //scaled problem in terms of the gcd of r_star - 1 and L
    mpz_t g;
    mpz_init( g );
    mpz_gcd( g, script_R, L );
    
    // the next line now has script_R holding R1 as described in Section 5.3.1 of ANTS 2024 paper
    mpz_divexact( script_R, script_R, g );
    
    // script_P = (P-1)/g
    mpz_t script_P;
    mpz_init_set( script_P, P);
    mpz_sub_ui( script_P, script_P, 1);
    mpz_divexact( script_P, script_P, g);
    
    // script_L = L/g
    mpz_t script_L;
    mpz_init_set( script_L, L);
    mpz_divexact( script_L, script_L, g);
    
    mpz_t div_bound;
    mpz_init_set( div_bound, script_P );
    mpz_sqrt( div_bound, div_bound );
    
    mpz_t divisor;
    mpz_set( divisor, script_R );
    
    // This constructs the smaller divisor and we should prevent P*R from exceeding B
    while( mpz_cmp( divisor, div_bound) <= 0 )
    {
        if( mpz_divisible_p( script_P, divisor ) )
        {
            mpz_mul( R, divisor, g);
            mpz_add_ui( R, R, 1);
            if( mpz_probab_prime_p( R, 0 ) != 0 )
            {
                // we might do a bounds check (?)
                // test that P*R is a CN
            }
        }
        mpz_add( divisor, divisor, script_L );
    }
    
    // Now we construct hte larger divisor and we should also prevent P*R from exceeding B
    // E.g, it should be impossible to enter this section if P > B^(2/3)
    
    // set R to be R_2 in section 5.3.2
    mpz_invert( script_R, script_R, script_L);
    mpz_mul( script_R, script_P, script_R );
    mpz_mod( divisor, script_R, script_L );
       
    while( mpz_cmp( divisor, div_bound) <= 0 )
    {
        if( mpz_divisible_p( script_P, divisor ) )
        {
            mpz_divexact( R, script_P, divisor);
            mpz_mul( R, R, g);
            mpz_add_ui( R, R, 1);
            if( mpz_probab_prime_p( R, 0 ) != 0 )
            {
                // we might do a bounds check (?)
                // test that P*R is a CN
            }
        }
        mpz_add( divisor, divisor, script_L );
    }
    
    mpz_clear( divisor );
    mpz_clear( div_bound );
    mpz_clear( script_R );
    mpz_clear( script_L );
    mpz_clear( script_P );
    mpz_clear( g );
    mpz_clear( R );

}


// while R is passed as mpz_t type
// the assumption is that R < 2^64
// NEED TO DO - incorporate append bound, see comments below
bool Preproduct::CN_factorization( mpz_t& n, mpz_t& R )
{
    std::queue<uint64_t> R_composite_factors;
    std::vector<uint64_t> R_prime_factors;
    
    uint64_t r64_factor =  mpz_get_ui( R );
    mpz_probab_prime_p( R, 0 ) == 0 ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
    
    mpz_t fermat_result;
    mpz_t odd_part;
    mpz_t base;
    mpz_t gcd_result;
    
    mpz_init( gcd_result );
    mpz_init_set_ui( base, 2);
    
    mpz_init( odd_part );
    mpz_init_set( fermat_result, n );
    // now using fermat_test to hold n-1
    mpz_sub_ui( fermat_result, fermat_result, 1 );
    uint16_t pow_of_2 = (uint16_t) mpz_scan1( fermat_result, 0 );
    mpz_fdiv_q_2exp( odd_part, fermat_result, pow_of_2 );
    
    
    
    uint16_t start_size;            // counter to empty queue
    uint16_t i = 0;                 // counter for powers of 2
    bool is_fermat_psp = true;      // initialized as true b/c n is a base-2 fermat psp

    // this while loops iterates over Fermat bases, starts with 2 and then odd integers
    mpz_powm( fermat_result, base, odd_part, n);
    mpz_set_ui( base, 1);
    while( !R_composite_factors.empty() && is_fermat_psp )
    {
        i = 0;    //counter for powers of 2
        // this while loop iterates over algebraic factors (b^(d * 2^i) + 1)
        // early abort?  no need to check all factors
        while( !R_composite_factors.empty() && i < pow_of_2 )
        {
            start_size = R_composite_factors.size();
            // this for loop iterates over composite factors
            for( int j = 0; j < start_size; j++ )
            {
                r64_factor = R_composite_factors.front();
                R_composite_factors.pop();
                mpz_set_ui( R, r64_factor );
                mpz_add_ui( gcd_result, fermat_result, 1); // gcd_result hold b^(d*2^i) + 1
                mpz_gcd( gcd_result, gcd_result, R);

                if( mpz_cmp(gcd_result, R) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
                {
                    // this gcd split r_factor.  Letting g = gcd( gr, rf), then this splits rf into, g and rf/g.
                    // could break here if rf or rf/g < append_bound (exiting early if not)
                    r64_factor = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0 ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
                    mpz_divexact(gcd_result, R, gcd_result );
                    r64_factor = mpz_get_ui( gcd_result );
                    mpz_probab_prime_p( gcd_result, 0 ) == 0  ? R_composite_factors.push( r64_factor ) : R_prime_factors.push_back( r64_factor );
                }
                else // r_factor was not factored, so put it back in the queue
                {
                    R_composite_factors.push( r64_factor );
                }
            }
            i++;
            mpz_powm_ui( fermat_result, fermat_result, 2, n );
        }

        // get new base
        mpz_add_ui( base, base, 2);
        mpz_powm( fermat_result, base, n, n);
        is_fermat_psp = ( mpz_cmp( fermat_result, base ) == 0 );
        mpz_powm( fermat_result, base, odd_part, n);
    }

    if( R_composite_factors.empty( ) && is_fermat_psp )
    {
        // could compare the least element in the sort
        std::sort ( R_prime_factors.begin(), R_prime_factors.end() );
        
        // should be in the below - removing to get more output
        // R_prime_factors[0] > append_bound &&
        if( appending_is_CN( R_prime_factors ) )
        {
            std::cout<< "          THIS IS A CARMICHAEL NUMBER     " << std::endl;
            gmp_printf ("%Zd = ", n );
            // use iterators !
            for( int i = 0; i < P_primes.size() ; i++ )
                std::cout << " " << P_primes[ i ] ;
            for( int i = 0; i < R_prime_factors.size(); i ++ )
                std::cout << " " << R_prime_factors[i];
            std::cout << std::endl;
        }
    }
        
    mpz_clear( fermat_result );
    mpz_clear( gcd_result );
    mpz_clear( base );
    mpz_clear( odd_part );
    
    return is_fermat_psp;
}



int main(void) {
    
    Preproduct P0;
    // P0.initializing( 599266767, 890750, 991 );
    P0.initializing( 6682828353, 2289560, 13 );
    
    std::cout << "LP for this is " << sizeof( unsigned long int ) << std::endl;
    
    std::cout << "Initializing P : " ;
    gmp_printf ("%Zd = ", P0.P );
    
    for( int i = 0; i < P0.P_primes.size(); i++)
    {
        std::cout << P0.P_primes[i] << " * "  ;       
    }
    std::cout << P0.P_primes[P0.P_primes.size() - 1 ] << std::endl ;
    
    std::cout << "Initializing Lambda : " ;
    gmp_printf ("%Zd = ", P0.L );
    
    for( int i = 0; i < ( P0.L_primes.size() - 1 ); i++)
    {
        std::cout << P0.L_primes[i] << " * "  ;
    }
    std::cout << P0.L_primes[P0.L_primes.size()  - 1 ] <<  std::endl ;

   
    /*  A good test case
    Preproduct P1;
    P1.initializing( 515410417841, 115920, 50 );
    P1.CN_search();
    */
     
    // the triple 132582235, 91872, 941 is in the working jobs with the new rule
    // this means we have between 2 and 6 primes to append to find a CN
    // for all primes in ( 941, (B/P)^(1/3) ) = (941, 196112) we append a prime and do CN_search
    // for all primes in ( 196112, (B/P)^(1/2) ) = ( 196112, 86847502) we append a prime and look for exactly one prime
    // this represents the first two CN_xearches invoked by the above:
    
    
    uint64_t P = 132582235;
    uint64_t L = 91872;
    
    uint64_t P1;
    uint64_t L1;
    
    uint64_t p64 = 947;
    mpz_t p ;
    mpz_init_set_ui( p, p64 );
    
    Preproduct PPP;
    PPP.initializing( 132582235, 91872, 900  );
   
    // this example is *not* ideal, we refactor P by trial division every time
    // we need to have p2-1 in factored form in order to use appending call
    // as you can see below, for ease, I'm just using gmp's next prime
    // and am just re-initializing each time.
    Preproduct P_testing;

    while( p64 < 196112 )
    {
        mpz_nextprime( p, p);
        p64 = mpz_get_ui( p );
        if( PPP.is_admissible_modchecks( p64 ) )
        {
            //std::cout << "I found an admissible prime " << p64 << std::endl;
            P1 = P*p64;
            L1 = L*((p64-1)/std::gcd( L , p64-1 ) );
            P_testing.initializing ( P1, L1, p64 );
            P_testing.CN_search();
            
        }
    }
    
    
    

    /*
     This run completes in 4 minutes and 30 seconds on thomas.butler.edu and finds these:
     433493717815335774335905 =  5 17 23 73 929 2377 7129 10267 18793
     72425097332690148535105 =  5 17 23 73 929 5347 577 673 263089
     433493717815335774335905 =  5 17 23 73 929 7129 2377 10267 18793
     433493717815335774335905 =  5 17 23 73 929 10267 2377 7129 18793
     433493717815335774335905 =  5 17 23 73 929 18793 2377 7129 10267
     725906640592907462305 =  5 17 23 73 929 21577 643 394633
     
     The fact that these are found duplicated only indicates our failure to implement the append_bound-related checks in CN_factorization
     */
   
    
    
    mpz_clear( p );
    
    // P0.CN_search(1873371784);
    //P0.CN_search(149637241475922);
    
    return 0;
}

