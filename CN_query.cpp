// Answers the query "Is n a Carmichael number?" in 
// in deterministic polynomial time.
// See Algorithm 1 and Theorem 9 of draft
// Relies on GMP's B-PSW "primality" test
// A "no" answer most likely returns a Fermat b that witnessed in is composite
// A "yes" answer returns the complete prime factorization for verification with Korselt's
//        A "no" answer can be found with Korselt's criterion, too.  

#include <cstdint>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <algorithm>
#include <vector>

int main()
{
    std::string string_number;
    std::cout << "This relies on GMP's B-PSW primality test and will answer the query:" << std::endl;
    std::cout << "Is n a Carmichael number?" << std::endl;
    std::cout << "Enter n: ";
    std::cin >> string_number; // Reads only the first word
	        
	mpz_class n;
	n = string_number;
	
    if( mpz_probab_prime_p( n.get_mpz_t(), 10 ) > 0 )
    {
        std::cout << "GMP detects that n is a (probable) prime." << std::endl;
    }
    else
    {
        std::vector< mpz_class > current_composites;
        std::vector< mpz_class > next_composites;
        std::vector< mpz_class > primes;
        std::vector< mpz_class > fermat_base_powers;
		        
        current_composites.push_back( n );
        
        mpz_class odd_part;
        odd_part = n - 1;
        // count powers of 2 dividing n-1:
        uint64_t pow_of_2 = (uint64_t) mpz_scan1( odd_part.get_mpz_t() , 0);
        // remove these powers of 2 to get the odd part:
        mpz_cdiv_q_2exp( odd_part.get_mpz_t(), n.get_mpz_t(), pow_of_2 );
        
        bool is_fermat_psp = true;
        
        mpz_class fermat_base = 1;
        
        while( !current_composites.empty() && is_fermat_psp )
        {
        	mpz_nextprime( fermat_base.get_mpz_t(), fermat_base.get_mpz_t() );
			mpz_class fermat_result;
			// compute fr = b^(odd_part) mod n
			mpz_powm( fermat_result.get_mpz_t(), fermat_base.get_mpz_t(), odd_part.get_mpz_t(), n.get_mpz_t() );
			fermat_base_powers.clear();
			for( uint64_t i = 0; i < pow_of_2; i++ )
			{
				fermat_base_powers.push_back( fermat_result );
				// compute fr = fr^2 mod n, this is b^(2^i * odd_part ) mod n
				mpz_powm( fermat_result.get_mpz_t(), fermat_result.get_mpz_t(), 2, n.get_mpz_t() );
			}
			// fr now holds b^(n-1) mod n.
			fermat_result = fermat_result * fermat_base % n;
			is_fermat_psp = ( fermat_result	== 	fermat_base );

			if( is_fermat_psp )
			{
				for( auto m : current_composites )
				{
					for( auto X : fermat_base_powers )
					{
						mpz_class g, temp;
						temp = X + 1;
						mpz_gcd( g.get_mpz_t(), m.get_mpz_t, temp.get_mpz_t() ); 
						m = m/g;
						if( g > 1 )
						{
							 mpz_probab_prime_p( g.get_mpz_t(), 10 ) == 0 ? next_composites.push_back( g ) : primes.push_back( g );
						}
					}
					if( m > 1 )
					{
						 mpz_probab_prime_p( m.get_mpz_t(), 10 ) == 0 ? next_composites.push_back( m ) : primes.push_back( m );
					}
				}
			}
			// we have completed a round of fermat testing and factoring
			current_composites.clear();
			std::swap( next_composites, current_composites );
        }
		// We have exited the while loop.  So is_fermat_psp is false or current_composites is empty
		if( !is_fermat_psp )
		{
			 std::cout << "the Fermat test detected your number is composite (and so, not a Carmichael number) with the base" << fermat_base << std::endl;
		}
		else
		{
			// Invoke Korselt now:
			mpz_class L = 1;
			mpz_class L_of_p;
			std::cout << n << " = 1 " ;
			for( auto p : primes )
			{
				std::cout << " * "  << p ;
				L_of_p = p - 1;
				mpz_lcm( L.get_mpz_t(), L.get_mpz_t(), L_of_p.get_mpz_t() )
			}
			std::endl;
			std::cout << "lambda( " << n << " ) = " << L << std::endl;
			L_of_p = ( n - 1 ) % L;
			std::cout << " (n-1) mod lambda(n) = " << L_of_p << std::endl;
			
			if( 0 == L_of_p )
			{
				std::cout << "n is a Carmichael number" << std::endl;
			}
			else
			{
				std::cout << "n is not a Carmichael number" << std::endl;
			}
				
		}
		
    }
}

