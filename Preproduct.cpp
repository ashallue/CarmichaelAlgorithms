#include "Preproduct.h"
#include <algorithm>
#include <iostream>
#include <queue>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>
#include <cstddef>
#include <boost/dynamic_bitset.hpp>

static_assert(sizeof(unsigned long) == 8, "unsigned long must be 8 bytes.  needed for mpz's unsigned longs to take 64 bit inputs in various calls.  LP64 model needed ");

// redo these if necessary

// hard code sqrt( B )
// we have choosen B = 10^24
#define SQRT_BOUND 1'000'000'000'000
// largest prime <= sqrt( B / X ) = 10^8
// because X = 10^8
#define DEFAULT_MAX_PRIME_BOUND 100'000'000
// there are 5761455 primes less than 10^8
#define PRIME_COUNT 5761455

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
    // So, we have an arry in sorted order and we need to insert a new element
    // copy the info from PP into the new product until we find where the new info goes
    // put the new info in
    // resume copying the info from PP if needed
    else
    {
        int i = 0;
        while( PP.next_inadmissible[i] < temp_next_inad && i < PP.len_appended_primes )
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


void Preproduct::CN_search(  )
{
    const uint32_t cache_bound = 150'000;
    
    mpz_t bound;
    mpz_init_set_ui( bound, 10 );
    mpz_pow_ui( bound, bound, 24 );
     
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
    
    mpz_t mpz_prime;
    mpz_init( mpz_prime );
    
    mpz_t small_prime;
    mpz_init( small_prime );
    
    mpz_t lifted_L;
    mpz_init_set( lifted_L, L );
    
    // position i has the truth value of the statement "(2i + 3) is prime"
    std::bitset<256> small_primes{"0010010100010100010000010110100010010100110000110000100010100010010100100100010010110000100100000010110100000010000110100110010010010000110010110100000100000110110000110010010100100110000110010100000010110110100010010100110100110010010110100110010110110111"};
    // fix append_bound to be consistent with the bitset
    uint64_t sieve_index_bound = std::min( (append_bound - 3)/2, (uint64_t) 512);
    
    // remove primes dividing L from the bitset
    uint16_t i = 1;
    while( i < L_len &&  L_distinct_primes[i] < 512 )
    {
        small_primes [ ( L_distinct_primes[i] - 3) /2 ] = 0;
        i++;
    }
    
    mpz_t cmp_bound;
    mpz_init( cmp_bound );
    mpz_cdiv_q( cmp_bound, bound, PL);
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
            mpz_cdiv_q( cmp_bound, bound, PL);
        }
        prime_index++;
    }

    uint32_t cmp_bound32 = mpz_get_ui( cmp_bound );
    boost::dynamic_bitset<> spoke_sieve( cmp_bound32 );
    
    for( auto p : primes_lifting_L )
    {
        std::cout << p << " " ;
    }
    
    // This is the wheel and a creates r_star + m*L
    // The wheel is r_star + m*L
    for( uint64_t m = 0; m < L_lift; m ++)
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
            while( mpz_cmp( n , bound ) < 0 )
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
                            // Call Fermat factorization here.
                            gmp_printf( "n = %Zd = %Zd * %Zd is a base-2 and base-3 Fermat psp. \n", n, P, R);
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
    mpz_clear( bound );
    mpz_clear( lifted_L );
        
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

std::vector< primes_stuff > Preproduct::primes_admissible_to_P( )
{
    std::vector< primes_stuff > return_vector;
    
    // there are 5761455 primes less than 10^8
    // we set aside the appropriate space
    // do this or let the space be dynamically allocated based on the computed prime_bound?
    return_vector.reserve( PRIME_COUNT );
    
    // a different way to do the below would be to
    // test if P > 10^8 first
    // only when P > 10^8 would prime_bound have a value less than
    mpz_t sqrt_P;
    mpz_init( sqrt_P );
    mpz_sqrt( sqrt_P, P );
    
    int64_t prime_bound;
    mpz_export( &prime_bound, 0, 1, sizeof(uint64_t), 0, 0, sqrt_P);
    prime_bound = SQRT_BOUND / prime_bound;
    prime_bound = std::min( prime_bound, (int64_t) DEFAULT_MAX_PRIME_BOUND );
    
    /*
     initialize some kind of factor sieve up to prime_bound
     use primes dividing P for admissibility checks
     e.g. for q in the factor sieve make sure 1 != q mod p for each p dividing P
     write each admissible prime into the prime_stuff struct formate
     push_back said prime_stuff
     */
    
    mpz_clear( sqrt_P );
    return return_vector;
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

/* Factor a Fermat pseudoprime n.  Fermat check not performed, just assumed.
   Prime, composite factors placed into appropriate vectors.
*/
void Preproduct::fermat_factor(uint64_t n, std::queue<uint64_t>& comp_factors, std::vector<uint64_t>& prime_factors, mpz_t& strong_result)
{
    // we are factoring n, so push n onto the composite queue, and clear the prime factors list
    comp_factors.push( n );
    prime_factors.clear();

    // temp variable to hold the factor pulled off the queue.  Then will be converted to mps_t n_factor.
    uint64_t temp;
    mpz_t n_factor;
    mpz_init( n_factor );

    // storage for the gcd result
    mpz_t gcd_result;
    mpz_init( gcd_result );
    
    int start_size = comp_factors.size();
    // use a for loop to go through all factors that are currently in the queue
    for( int j = 0; j < start_size; j++ )
    {
        // get element out of queue and put into mpz_t
        // first time through, this is just n
        temp = comp_factors.front();
        comp_factors.pop();
        
        mpz_set_ui( n_factor, temp );
          
        // check gcd before prime testing
        // strong_result holds the algebraic factor assoicated with b^((n-1)/2^e) + 1
        mpz_add_ui( strong_result, strong_result, 1);
        mpz_gcd( gcd_result, strong_result, n_factor);
        
        // check that gcd_result has a nontrivial divisor of n_factor
        // could probably be a check on result1 = +/- 1 mod n
        // before computing the gcd
        if( mpz_cmp(gcd_result, n_factor) < 0 && mpz_cmp_ui(gcd_result, 1) > 0 )
        {
          // will need to add a check about a lower bound on these divisors
          mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result);
          ( mpz_probab_prime_p( gcd_result, 0 ) == 0 ) ? comp_factors.push( temp ) : prime_factors.push_back( temp );
          mpz_divexact(gcd_result, n_factor, gcd_result );
          mpz_export( &temp, 0, 1, sizeof(uint64_t), 0, 0, gcd_result);
          ( mpz_probab_prime_p( gcd_result, 0 ) == 0 ) ? comp_factors.push( temp ) : prime_factors.push_back( temp );
        }
        else // n_factor was not factored, so it is prime or composite
        {
          ( mpz_probab_prime_p( n_factor, 0 ) == 0 ) ? comp_factors.push( temp ) : prime_factors.push_back( temp );
        }
    }
    // if R_composite is empty, check n is CN *here*
    // output lines below are temporary and meant for debugging
    gmp_printf ("factoring n = %Zd", n);
    std::cout << " has " << comp_factors.size() << " composite factors and " << prime_factors.size() << " prime factors." << std::endl;

    // need to clear mpz vars
}

/* Check whether n is a Fermat pseudoprime to the base b.  Returns bool with this result.
   Additionally, sets strong_result variable to b^((n-1)/2^e) + 1
   Notes this function returns true for prime n.
 
*/
bool Preproduct::fermat_test(mpz_t& n, mpz_t& b, mpz_t& strong_result)
{
     
    // create a variable for n-1, then compute the largest power of 2 that divides n-1
    mpz_t nminus;
    mpz_init( nminus );
    mpz_sub_ui( nminus, n, 1);
    // counts the number of 0's that terminate in nminus, i.e. e such that 2^e || n-1
    uint32_t exp_on_2 = (uint16_t) mpz_scan1( nminus , 0);
    uint32_t pow_of_2 = ( 1 << exp_on_2 );
    // QUESTION: should pow_of_2 be larger than 32 bits?
    // ANSWER:  no.  Indeed, a char would work.  everything in this computation can fit in an uint128_t
    // which means that exponent is bounded by 128 which is a 7 bit number.
    
    
    // stores the exponent we will apply to the base
    mpz_t strong_exp;
    mpz_init( strong_exp );

    // stores the result b^(n-1) mod n, i.e. the Fermat exponent
    mpz_t fermat_result;
    mpz_init( fermat_result );
    
    // set up strong base:  truncated divsion by 2^e means the exponent holds (n-1)/(2^e)
    mpz_tdiv_q_2exp( strong_exp, nminus, exp_on_2 );
    // we use prime divisors of L as the Fermat bases
    mpz_powm( strong_result,  b,  strong_exp, n); // b^( (n-1)/(2^e) ) mod n
    mpz_powm_ui( fermat_result,  strong_result, pow_of_2, n); // b^( (n-1)/(2^e)) )^(2^e) = b^(n-1)

    //std::cout << "strong_exp = " << mpz_get_str(NULL, 10, strong_exp) << " strong_result = " << mpz_get_str(NULL, 10, strong_result);
    //std::cout << " fermat_result = " << mpz_get_str(NULL, 10, fermat_result) << "\n";

    bool is_psp = mpz_cmp_ui( fermat_result, 1 ) == 0;

    // deallocating mpz vars created in this function
    mpz_clear( nminus );
    mpz_clear( strong_exp );
    mpz_clear( fermat_result );
    
    return is_psp;   
}


int main(void) {
    
    Preproduct P0;
    // P0.initializing( 599266767, 890750, 991 );
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

    std::cout << "Testing is_fermat\n";
    mpz_t n;
    mpz_init(n);
    mpz_t base;
    mpz_init(base);
    mpz_set_ui(base, 2);
    mpz_t strong_result;
    mpz_init(strong_result);
    
    for(int i = 100; i < 10000; i++){
        mpz_set_ui(n, i);
        bool is_psp = P0.fermat_test(n, base, strong_result);
        if(is_psp && mpz_probab_prime_p( n, 0 ) == 0 ){
            std::cout << i << " is a pseudoprime\n";
        }
    }
    
    
    // P0.CN_search(1873371784);
    //P0.CN_search(149637241475922);
    
    return 0;
}

