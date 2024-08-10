#include <gmpxx.h>
#include <iostream>

int main()
{
	mpz_class n, nm1,  offset, output, base, result;
	base = 2;
	// offset is  P*lambda(P)
	// n is P*r where r is the least positive residue of P^{-1} mod lambda(P)
	// I used the online magma calculator to compute these and pasted them in
	// http://magma.maths.usyd.edu.au/calc/
	n = "21415192274146561";
	offset = "143890337981477760";
	// n = "370031963673555618862708354823205121";
	// offset = "1329223147719610607895840242550071040";
	// n+= 2'000'000'000*offset;
	nm1 = n - 1;


	for( int64_t i = 0; i < 1'000'000'000; i++)
	{
		mpz_powm( result.get_mpz_t(),  base.get_mpz_t(),  nm1.get_mpz_t(), n.get_mpz_t());
		if( cmp ( result, 1 ) == 0 )
		{
			std::cout << " i = "  << i << " and n = " << n << std::endl;
		}
		n+= offset;
		nm1 += offset;
	}
}
