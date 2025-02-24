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

    // constructor always sets to 10^24
    mpz_t BOUND;
  
    // stores P and the primes dividing p in sorted order
    mpz_t P;
    std::vector< uint64_t > P_primes;
    
    // previous versions had extra data associated to P_primes
    // would store candidate for inadmissibility so that admissibiliy checks were inequality
    // if the inequality failed, an "easy" update would restore the data strucutre updating the least candidate
    // this is still viable - a hybrid of mod checks for small primes and the data structure for larger primes
   
    // stores L
    mpz_t L;
    
    // previous versions had extra data assocaited with L
    // if the complete prime factorization of L and p-1 are known
    // then LCM( L, p - 1) may be computed via a merge
    // and the euclidaean algorithm need not be invoked
    // we might revisit this but in the present version we do not do this and removed this support

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
    void appending( Preproduct& PP, uint64_t prime );
    
    // member functions
       
    // done with no gcd check - checks if input is not 1 mod p for all p | P
    // consdier Lemire's fast remainders (?)
    // see also comment above
    bool is_admissible_modchecks( uint64_t prime_to_append );
    
    
    // uses a rule to decide whether to call CN_search or
    // initialize an incremental sieve and do prime-by-prime appending
    // for primes q in (append_bound, (B/P)^(1/3) ) we recursively call this on preproduct Pq
    // for primes q in ((B/P)^(1/3) , (B/P)^(1/2) ) we do the factorization on Pq-1
    void complete_tabulation( std::string cars_file );
    
    
    // This does *not* create a Preproduct structure
    bool appending_is_CN( std::vector< uint64_t >&  primes_to_append );

    // finds all R = ( P^{-1} mod L ) + k*L satisfying P*R < B
    // does a sieve on k values and removes candidates where
    //      R would be divisibly be a prime less than or equal append_bound
    //      or R would be divisible by a prime inadmissible to P
    // checks that each candidate is a Fermat psuedoprime
    // calls CN_factorization when P*R ius a base 2 and base 3 Fermat psp
    void CN_search( std::string cars_file  );
    
    
    // For some preproducts, we will still only need to find exactly one prime
    // Implements the modified trial division of Serction 5.3 ANTS 2024 paper
    // We could consider two version:
    // Version 1 - a custom version of CN_search that only outputs primes (consider this the bounded factorization)
    // Version 2 - the Lenstra Divisors in residue class algorithm
    // Contrary-wise:
    // Our implementation accounts for the bounded case (although does not do prime sieving when bounded)
    // and uses residue class information in a straighforward way (easier to check the correctness)
    void completing_with_exactly_one_prime( );
    
    // check that L exactly divides P - 1
    // we might not need this
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
    bool CN_factorization( mpz_t& n, mpz_t& R, std::vector<uint64_t>& R_prime_factors, std::string cars_file );


};

#endif
