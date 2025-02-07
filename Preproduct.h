#ifndef PREPRODUCT_H
#define PREPRODUCT_H

#include <algorithm>
#include <vector>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>
#include <queue>
#include <vector>

class Preproduct{
    

public:

    mpz_t BOUND;  // constructor always sets to 10^24
  
    mpz_t P;
    std::vector< uint64_t > P_primes;
   
    // stores the distinct primes that divide L
    mpz_t L;
    std::vector< uint64_t > L_primes;

    // primes appended to P need to exceed this bound
    uint64_t append_bound;

    // constructor and destructor
    Preproduct();
    ~Preproduct();
    
    // initializing call
    // has to factor init_preproduct and init_LofP
    void initializing( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound );
    
    // appending call
    // assume we have an admissible prime to append.
    // contains a merge computation of LCM( lambda(PP), p-1 )	
    void appending( Preproduct* PP, uint64_t prime, std::vector< uint64_t >& distinct_primes_dividing_pm1 );
    
    // member functions
       
    // done with no gcd check - checks if input is not 1 mod p for all p | P
    bool is_admissible_modchecks( uint64_t prime_to_append );
    
    
    // uses a rule to decide whether to call CN_search or
    // initialize an incremental sieve and do prime-by-prime appending
    // for primes q in (append_bound, (B/P)^(1/3) ) we recursively call this on preproduct Pq
    // for primes q in ((B/P)^(1/3) , (B/P)^(1/2) ) we do the factorization on Pq-1
    void complete_tabulation( );
    
    
    // This will compute L and P with gcd computations
    // does *not* create a Preproduct structure
    // future version should probably have a filestream argument
    // and have this method be void but write output to file
    bool appending_is_CN( std::vector< uint64_t >&  primes_to_append );
    
    // finds all R = ( P^{-1} mod L ) + k*L satisfying P*R < B
    // checks that each candidate is a Fermat psuedoprime
    // calls CN_factorization when P*R ius a base 2 and base 3 Fermat psp
    void CN_search( );

    
    // For some preproducts, we will still only need to find exactly one prime
    // B/(PL) is small, just call CN_search instead
    // one could work out the cross over depending on the implementation of the below
    // First attempt is the modification of trial division due to the residue class information
    // To do (?): Lenstra's divisor in residue class
    // To do : incorporate bounded search
    void completing_with_exactly_one_prime( );
    
    // check that L exactly divides P - 1
    // in the future modify to take filestream?
    bool is_CN( );

    
    // passes n and R where n = P*R
    // uses strong Fermat primality tests for fast factorization
    // if a simple Fermat test detects a composite number, n cannot be a CN
    // otherwise, it is likely that n will be detected as composite by the strong test
    // in which case, we may use the strong test to split R
    // see:  https://crypto.stackexchange.com/questions/5279/carmichael-number-factoring
    // in the above link they claim a "typical" number is detected as compoiste with probaility  > 3/4
    // but that for CN the probability is 7/8
    // I do not know where the 7/8 comes from
    bool CN_factorization( mpz_t& n, mpz_t& R, std::vector<uint64_t>& R_prime_factors);


};

#endif
