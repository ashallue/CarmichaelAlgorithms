#ifndef _ROLLSIEVE
#define _ROLLSIEVE

/*
Incrementally finds primes using the rolling sieve

Rollsieve r(start);     -- constructor, will enumerate primes > start
bool r.next();          -- true if next integer is prime; advances to next integer
uint64_t r.nextprime(); -- returns the next larger prime
uint64_t r.getn();      -- returns the current value of n

To list all primes from 1 to n, use this:

 Rollsieve r(1);
 for(uint64_t p=r.nextprime(); p<=n; p=r.nextprime()) { std::cout << p << std::endl; }

Credit: Jonathan Sorenson
*/


#include<cmath>
#include<cstdint>
#include<vector>


//using namespace std;

class Factorlist2
{
    public:
    static const int maxplen=15;
    uint64_t list;
    char plen;
    char ptr[maxplen];

    Factorlist2(): list(0), plen(0) {}

    inline void clear() { list=0; plen=0; }

    uint32_t get(uint32_t pos);

    uint32_t getbitlen(uint32_t pos);

    uint32_t bitlength(uint32_t x);
    void add(uint32_t p, uint32_t bitlen);
    
    inline void add(uint32_t p) { add(p,bitlength(p)); }
    inline void push(uint32_t p, int bitlen) { add(p,bitlen); }
    inline void push(uint32_t p) { add(p); }
    inline uint32_t pop() { return get( --plen ); }
    inline bool isempty() { return (plen==0); }
    inline uint32_t gettop() { return get(plen-1); }
    inline uint32_t gettopbitlen() { return getbitlen(plen-1); }
    inline void makeempty() { clear(); }

    void getlist(uint64_t n, std::vector<uint64_t> & plist);

};  // end of Factorlist2 class

class Rollsieve
{

    static const int primeslen = 168;
    static const uint16_t primesmax = 1000;
    static const uint16_t primes[168]; // initialized at the end of the file

    std::vector<Factorlist2> T;
    uint32_t delta;

    uint64_t n, s;
    uint32_t pos, r;

    public:
    Rollsieve(uint64_t start);

    inline uint64_t getn() { return n; }

    bool next();  // this code is nearly verbatim from the paper
  
    uint64_t nextprime() { while(!next()); return n-1; }
    
    bool isnextprime();

    void getlist(std::vector<uint64_t> & plist) { T[pos].getlist(n,plist); }
};


#endif
