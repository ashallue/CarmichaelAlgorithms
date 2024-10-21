#ifndef PREPRODUCT_H
#define PREPRODUCT_H


class Preproduct{
  public:

    const mpz_t P;
    const uint64_t P_primes[20];
    const uint16_t P_len;
    const uint64_t append_bound;  // primes appended to P need to exceed this bound

    // Information about L = CarmichaelLambda(P)
    const uint64_t L;
    const uint64_t L_distinct_primes[20];
    const uint16_t  L_exponents[20];   
    const uint16_t L_len;

    // 5 is chosen because that's the limit of the recursion depth.
    const uint16_t len_appended_primes;

    // these arrays are used to avoid gcd computations for admissibility checks
    // they are updated assuming primes are tested for admissibility in increasing order
    // These three arrays are the only non const type - requires updating
    uint64_t next_inadmissible[5]; 
    uint16_t mod_three_status[5];  
    uint64_t appended_primes[5];   


    // constructor and destructor here


    //
    bool is_addmissible( uint64_t prime_to_append );

};

#endif
