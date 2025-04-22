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
    
    // lower level improvements that are currently not implemented:
    // 1 - we can do this with uint128_t and avoid GMP entirely
    //      1a - the primary need is something faster than GMP to compute 2^(n-1) mod n
    //           this is the primary elimination check in CN_search
    //           for P = 515410417841, we did 16,737,417 modular exponentiations and only found 126 base 2 psps
    //      1b - while writing a new mdoular exponentiation routine, we could integrate CN_factorization into it
    //           I think (needs to be checked) that once we transform to the Montgomery representation, the gcds
    //           can be done in that representation.  So, we do not need to enter-exit as we do right now
    //      1c - most of the other number theoretic functions can be found in std or are more easily written
    // 2 - partitition the P_primes vector into by prime size
    //      2a - for small primes use fast remainders for admissibility, e.g.:
    //           https://arxiv.org/abs/1902.01961
    //      2b - for other primes, keep store the next least inadmissible candidate
    //           then admissbility is checking an inequality (updating data-structure if not)
    //           by chooseing the not-small primes large enough the data structure should not be updated often
    // 3 - store explicit prime factorization of L and q-1, both primes and their exponents
    //           updating L can then be done with a merge of
    //           the merge and multiplications have to be cheaper than the Euclican algorithm approach with LCM( L ,q-1)
    // 4 - At all stages we compute P^{-1} mod L.  Is it possible to incrementally update this quantity?
    //           That is, given P, P^{-1} mod L, q a prime admissible to P (and whatever prime factorizations we might want)
    //           can we compute (Pq)^{-1} mod LCM( L , q-1 ) faster than the obvious method?
    // 5 - A better implementation (or two) of completing_with_exactly_one_prime
    //           5a - For case 3 of complete_tabulation, the expection is that this results in the bounded factorization
    //           5b - For case 1 of complete_tabulation, the expecation is tha this results in the unbounded factorization
    //                  so, maybe this can be done faster by using Lenstra's version rather than our modified trial divison

    // constructor always sets to 10^24
    mpz_t BOUND;
  
    // stores P and the primes dividing p in sorted order
    mpz_t P;
    std::vector< uint64_t > P_primes;
      
    // stores L
    mpz_t L;
    
    // primes appended to P need to exceed this bound
    uint64_t append_bound;

    // constructor and destructor
    Preproduct();
    ~Preproduct();
    
    // initializing call
    // has to factor init_preproduct
    void initializing( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound );
    
    // appending call
    // we assume the prime is admissible to PP
    void appending( Preproduct& PP, uint64_t prime );
    
    // done with no gcd check - checks if input is not 1 mod p for all p | P
    bool is_admissible_modchecks( uint64_t prime_to_append );

    void CN_multiples_of_P( std::string cars_file );
              
    // This does *not* create a Preproduct structure
    bool appending_is_CN( std::vector< uint64_t >&  primes_to_append, std::string cars_file  );

    // finds all R = ( P^{-1} mod L ) + k*L satisfying P*R < B
    // does a sieve on k values and removes candidates where
    //      R would be divisibly be a prime less than or equal append_bound
    //      or R would be divisible by a prime inadmissible to P
    // checks that each candidate is a Fermat psuedoprime
    // calls CN_factorization when P*R ius a base 2 and base 3 Fermat psp
    void CN_search( std::string cars_file  );
    
    // For some preproducts, we will still only need to find exactly one prime
    // Implements the modified trial division of Serction 5.3 ANTS 2024 paper
    void completing_with_exactly_one_prime( std::string cars_file );
    
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
    bool CN_factorization( mpz_t& n, mpz_t& R, std::string cars_file );

};

#endif
