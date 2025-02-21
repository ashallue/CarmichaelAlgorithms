#include "rollsieve.h"
#include <iostream>
#include <algorithm>
#include <gmp.h>
#include <numeric>

int main(void) {
    
    std::vector< uint64_t > pm1_factors;

    Rollsieve r(1);

    uint64_t p = r.getn();
    
    while( p < 400 )
    {
        r.next();
        p = r.getn();
    }
    
    
    while( p < 600 )
    {
        r.getlist( pm1_factors );
        std::cout << "prime factors of " << p << " are " ;
        std::sort( pm1_factors.begin(), pm1_factors.end() );
        for( auto f : pm1_factors )
        {
            std::cout << f << " " ;
        }
        std::cout << std::endl;

        r.next();
        p = r.getn();
    }
    


    return 0;
}


