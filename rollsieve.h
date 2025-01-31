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

*/

//#include "int.h"
#include<cmath>
#include<cstdint>
#include<vector>
//#include"factorlist2.h"

using namespace std;

class Rollsieve
{

    class Factorlist2
    {
        public:
        static const int maxplen=15;
        uint64_t list;
        char plen;
        char ptr[maxplen];

        Factorlist2(): list(0), plen(0) {}

        inline void clear() { list=0; plen=0; }

        inline uint32_t get(uint32_t pos)
        { // assumes 0<=pos<plen
            uint32_t left;
            left= (pos>0)?ptr[pos-1]:0;
            uint64_t copy;
            copy=(list >> left);
            return (uint32_t)
            ((copy & ( (1ul<<(ptr[pos]-left) ) -1 ))<<1)|1ul;
        }

        inline uint32_t getbitlen(uint32_t pos)
        { // assumes 0<=pos<plen
            uint32_t left;
            left= (pos>0)?ptr[pos-1]:0;
            return ptr[pos]-left+1;
        }

        static uint32_t bitlength(uint32_t x)
        {
            uint32_t len=0;
            while((1ul<<len) < x) len++;
            return len;
        }

        inline void add(uint32_t p, uint32_t bitlen)
        {
            uint32_t pos=plen;
            plen++;
            uint32_t left= (pos>0)?ptr[pos-1]:0;
            uint64_t pbits=(p>>1);
            list = list | (pbits<<left);
            ptr[pos]=left+bitlen-1;
        }
        
        inline void add(uint32_t p) { add(p,bitlength(p)); }
        inline void push(uint32_t p, int bitlen) { add(p,bitlen); }
        inline void push(uint32_t p) { add(p); }
        inline uint32_t pop() { return get( --plen ); }
        inline bool isempty() { return (plen==0); }
        inline uint32_t gettop() { return get(plen-1); }
        inline uint32_t gettopbitlen() { return getbitlen(plen-1); }
        inline void makeempty() { clear(); }

        void getlist(uint64_t n, vector<uint64_t> & plist)
        {
            plist.clear();

            // is n even? 2 is not in list
            if((n&1ul)==0)
            {
                plist.push_back(2); n>>=1;
                while((n&1ul)==0) n>>=1;
            }

            uint32_t p;
            for(int i=0; i<plen; i++)
            {
                p=get(i);
                plist.push_back(p);
                while(n%p==0) n=n/p;
            }
            if(n>1) plist.push_back(n);
        }

    };  // end of Factorlist2 class

    static const int primeslen = 168;
    static const uint16_t primesmax = 1000;
    static const uint16_t primes[168]; // initialized at the end of the file

    vector<Factorlist2> T;
    uint32_t delta;

    uint64_t n, s;
    uint32_t pos, r;

    public:
    Rollsieve(uint64_t start)
    {
        if(start<2) start=2;
        uint32_t sqrtstart = (uint64_t) sqrtl((long double) start);
        r=sqrtstart+1; s=r*r;
        delta=r+2;
        pos=0; n=start;
        
        T.reserve(2*delta+1000); // some wiggle room
        T.resize(delta);
        for( uint32_t i = 0; i < delta; i++) T[i].makeempty();
        
        // small primes first - take from int.h array if possible
        // for(int i=0; i<primeslen && primes[i]<=sqrtstart; i++)
        for( uint32_t i=1; i<primeslen && primes[i]<=sqrtstart; i++)
        {
            uint32_t p=primes[i];
            uint32_t j=(p-(start%p))%p;
            T[j].push(p);
        }
        // now large primes
        if(sqrtstart>primesmax)
        {
            Rollsieve R(primesmax); // hooray for recursion
            uint32_t p=R.nextprime();
            while(p<=sqrtstart)
            {
                uint32_t j=(p-(start%p))%p;
                T[j].push(p);
                p=R.nextprime();
            }
        }
    }

    inline uint64_t getn() { return n; }

    bool next()  // this code is nearly verbatim from the paper
    {
        bool isprime=(n%2!=0 || n==2);  // we no longer store 2
        while(!T[pos].isempty())
        {
            uint32_t bitlen=T[pos].gettopbitlen();
            uint32_t p=T[pos].pop();
            T[(pos+p)%delta].push(p,bitlen);
            isprime=false;
        }
        T[pos].clear();
        if(n==s)
        {
            if(isprime)
            {
                isprime=false;
                T[(pos+r)%delta].push(r);
            }
            r++; s=r*r;
        }
        n++; pos=(pos+1)%delta;
        if(pos==0)
        {
            delta+=2;
            T.resize(delta);
            T[delta-1].makeempty();
            T[delta-2].makeempty();
        }
        return isprime;
    }
  
    uint64_t nextprime() { while(!next()); return n-1; }
    
    bool isnextprime()
    {
        // checks if n+1 is prime:
        // 1) n is even so that n+1 is odd
        // 2) T( pos + 1 ) needs to be empty
        // 3) n + 1 could be the square of the next unrecorded prime, check that it is not
        return ( (0 == n%2) && T[ (pos+1) % delta ].isempty() && ( s != ( n + 1 ) ) ) ;
    }

    void getlist(vector<uint64_t> & plist) { T[pos].getlist(n,plist); }
};

const uint16_t Rollsieve::primes[168] = {
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

#endif
