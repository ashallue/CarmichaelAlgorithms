#include "rollsieve.h"
#include <iostream>
#include <algorithm>
#include <gmp.h>
#include <numeric>

int main(void) {
    
    std::vector< uint64_t > pm1_factors;

    Rollsieve r(1);

    uint64_t p = r.getn();
    
    while( p < 100 )
    {
        std::cout << p << std::endl;
        if( r.isnextprime() )
        {
            std::cout << p + 1 << " is prime." << std::endl;
            r.getlist( pm1_factors );
            for( auto f : pm1_factors )
            {
                std::cout << f << " " ;
            }
            std::cout << std::endl;
        }
        r.next();
        p = r.getn();
    }
    


    return 0;
}


