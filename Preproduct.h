#ifndef PREPRODUCT_H
#define PREPRODUCT_H


class Preproduct{
  public:
    mpz_t P;
    uint64_t P_primes[20];
    uint8_t P_len;
    uint64_t append_bound;  // primes appended to P need to exceed this bound

    // Information about L = CarmichaelLambda(P)
    uint64_t L;
    uint64_t L_distinct_primes[20];
    uint8_t  L_exponents[20];   // no exponent will exceed 255
    uint8_t L_len;

    //  maybe we need a child class of this class or maybe it needs two such classes
    //  The first class is the initiating preproduct
    //  The other classes are the derived preproducts from appending primes
    //  attempting as a single class for the moment

    // 5 is chosen because that's the limit of the recursion depth.  
    uint8_t primes_appended;
    uint64_t next_inadmissible_prime[5];
    bool mod_three_status[5];
    uint64_t P_primes[5];

};

#endif
