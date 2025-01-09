#ifndef PREPRODUCT_H
#define PREPRODUCT_H

#include <algorithm>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>
#include <queue>
#include <vector>

// we could consider a re-write for L and prime_stuff
// we could only store the exponent for 2
// as it stands pm1_distinct_primes[ 0 ] alwaqys hold 2
// and L_distinct_primes[ 0 ] always holds 2
// we could generalize this to the other small primes
// but the guaranteed divisibility is only for 2

// http://www.s369624816.websitehome.co.uk/rgep/cartable.html
// the least CN with a fixed number of prime factors is known
// array lengths set at fixed length to accomodate this
// for B = 10^24 there is 1 CN with 14 prime factors
// could be set smaller since the last prime(s) to be
// appended are done so without forming explicitly forming a Preproduct object
#define MAX_PRIME_FACTORS 14

// the intended use of this is with a precomputation that limits
// the total number of primes that can be appended in a prime-by-prime way
// for this computation, we choose that bound to be 5
#define APPEND_LIMIT 5

// contains a prime p and the factorization information for p-1 = Lambda(p)
// prime_stuff is for the appended prime that come from
// sieving up to sqrt( B/P )
// a quick check of all cases for p < 10^8
// resulted in at most 8 distinct prime factors for p-1
// an alternative data-structure for prime_suff
// might be only store the exponent for small primes:  2 3 5 and 7
// and store all prime factors (duplicates included) for all other primes in a single array
// "merge" computation of lcm( L(P), p-1) would be a bit easier
#define L_PRIME_FACTORS 8

// consider a more compact storage structure
// see:  https://github.com/sorenson64/soespace/blob/main/soe.h
// and the paper:  https://arxiv.org/pdf/2406.09150
// we'd nee to expand to something like the below for each append call
// but the storage cost would be reduced by about 16 times
struct primes_stuff
{
    uint32_t prime;
    uint32_t pm1_distinct_primes[ L_PRIME_FACTORS ];
    uint16_t pm1_exponents[ L_PRIME_FACTORS ];
    uint16_t pm1_len;
};

class Preproduct{
    
	
public:

    mpz_t P;
    uint64_t P_primes[ MAX_PRIME_FACTORS ];
    uint16_t P_len;
    uint64_t append_bound;  // primes appended to P need to exceed this bound

    // Information about L = CarmichaelLambda(P)
    // should change L to mpz_t
    // uint64_t L;
    mpz_t L;
    uint64_t L_distinct_primes[ L_PRIME_FACTORS ];
    uint16_t L_exponents[ L_PRIME_FACTORS ];   
    uint16_t L_len;

    // two forms of initialization
    // 1) "intializing" preproduct from the precomputation phase
    // 2) "appending" to the initializing product when prime by prime is justified
    // in case 2, we count these appended primes
    uint16_t len_appended_primes;
    // these arrays are used to avoid gcd computations for admissibility checks
    // updated assuming primes are tested for admissibility in increasing order
    uint64_t next_inadmissible[ APPEND_LIMIT ];
    uint16_t mod_three_status[ APPEND_LIMIT ];  
    uint64_t appended_primes[ APPEND_LIMIT ];  

    // small primes that serve as Fermat bases
    static const uint32_t smallprimeslen = 168;
    static const uint32_t smallprimesmax = 1000;
    uint32_t smallprimes[smallprimeslen] = {
       2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
      31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
      73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
     127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
     179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
     233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
     283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
     353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
     419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
     467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
     547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
     607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
     661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
     739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
     811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
     877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
     947, 953, 967, 971, 977, 983, 991, 997
    };

    // constructor and destructor
    Preproduct();
    ~Preproduct();
    
    // initializing call
    // has to factor init_preproduct and init_LofP
    void initializing( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound );
    
    // appending call
    // assume we have an admissible prime to append.
    // contains a merge computation of LCM( lambda(PP), p-1 )	
    void appending( Preproduct PP, primes_stuff p );

    // member functions
    // done with no gcd check
    bool is_admissible( uint64_t prime_to_append );

    // This will compute L and P with gcd computations
    // does *not* create a Preproduct structure
    // future version should probably have a filestream argument
    // and have this method be void but write output to file
    bool appending_is_CN( std::vector< uint64_t >&  primes_to_append );
    
    // finds all R = ( P^{-1} mod L ) + k*L satisfying P*R < B
    // checks that each candidate is a Fermat psuedoprime
    // uses a stronger Fermat test to factor composite R
    // if R is fully factored (and has passed the Fermat tests)
    // then P*R is to be checked with Korselt's criterion  NYI
    // meant to be called when it is no longer efficient to do prime-by-prime appending 
    // this takes the bound on R as an argument which implies R <= (B/P) < 2^64
    // and that L < 2^64
    void CN_search( );

    
    // For some preproducts, we will still only need to find exactly one prime
    // B/(PL) is small, just call CN_search instead
    // one could work out the cross over depending on the implementation of the below
    // First attempt is the modification of trial division due to the residue class information
    // To do (?): Lenstra's divisor in residue class
    void completing_with_exactly_one_prime( );
    
    // finds all primes that are admissible to P
    // the intent is that this creates the vector that holds the primes
    // that are used with the appending method
    std::vector< primes_stuff > primes_admissible_to_P( );
    
    // check that L exactly divides P - 1
    // in the future modify to take filestream?
    // to return true, means we need to output which is the actual goal
    bool is_CN( );

    /* Factor n, storing the unique prime factors in the associated parameter.
       Technique is the Fermat method: if n passes the fermat test to a base b, the associated 
       strong test can be used to split.  If n fails a fermat test, quit and return false to signify not carmichael.
    */
    bool fermat_factor(uint64_t n, std::vector<uint64_t>& prime_factors);

    /* Check whether n is a Fermat pseudoprime to the base b.  Returns bool with this result.
       Additionally, sets strong_result variable to b^((n-1)/2^e)
       Note this function returns true for prime n.
    */
    bool fermat_test(uint64_t& n, mpz_t& b, mpz_t& strong_result);
};

#endif
