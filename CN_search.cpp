// compiled with  g++ CN_search.cpp -lgmpxx -lgmp -O3

#include <gmpxx.h>
#include <iostream>
#include <chrono>

int main()
{



	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	auto t1 = high_resolution_clock::now();

	// version 1
	// version 1 took 401k ms

	/*
	mpz_class P, L, r_star, n, nm1,  offset, output, base, result;

	// Hard coded for P = F_4 and L = Lambda(F_5) and r_star = P^{-1} mod L
	P = 4294967297;
	L = 33502080;
	r_star = 4986113;

	n = P*r_star;
	nm1 = n - 1;

	offset = P*L;
	base = 2;

	for( int64_t i = 0; i < 100'000'000; i++)
	{
		mpz_powm( result.get_mpz_t(),  base.get_mpz_t(),  nm1.get_mpz_t(), n.get_mpz_t());
		if( cmp( result, 1 ) == 0 )
		{
			std::cout << " i = "  << i << " and n = " << n << std::endl;
		}
		n+= offset;
		nm1 = n - 1;
	}
  */

	auto t2 = high_resolution_clock::now();
	duration<double, std::milli> ms_double = t2 - t1;
	std::cout << ms_double.count() << "ms\n";
	//version 2:  uses int128 where possible, converstion to mpz_t for modular exponentiation

	t1 = high_resolution_clock::now();

	__int128 P1 = 4294967297;
	__int128 L1 = 33502080;
	__int128 r_star1 = 4986113;
	__int128 n1 = P1*r_star1;
	__int128 offset1 = P1 *L1;
	__int128 n1m1 = n1 -1 ;

	mpz_t base1;
	mpz_t result1;
	mpz_t n_mpz;
	mpz_t nm1_mpz;

	mpz_init( base1 );
	mpz_init( result1 );
	mpz_init( n_mpz );
	mpz_init( nm1_mpz );

	mpz_set_ui( base1, 2);

	for( int64_t i = 0; i < 100'000'000; i++ )
	{
		mpz_import( n_mpz, 1, 1, 16, 0, 0, &n1 );
		mpz_import( nm1_mpz, 1, 1, 16, 0, 0, &n1m1 );

		mpz_powm( result1,  base1,  nm1_mpz, n_mpz);
		if( mpz_cmp_si( result1, 1 ) == 0 )
		{
			std::cout << " i = " <<  i << std::endl;
		}
		n1+= offset1;
		n1m1 = n1 - 1;
	}

	mpz_clear( base1 );
	mpz_clear( result1 );
	mpz_clear( n_mpz );
	mpz_clear( nm1_mpz );

	t2 = high_resolution_clock::now();
	ms_double = t2 - t1;
	std::cout << ms_double.count() << "ms\n";

	//version 3 only uses mpz_t

	t1 = high_resolution_clock::now();

	mpz_t P;
	mpz_init(P);
	mpz_set_ui( P, 4294967297);

	mpz_t L;
	mpz_init(L);
	mpz_set_ui( L, 33502080);

	mpz_t r_star;
	mpz_init(r_star);
	mpz_set_ui( r_star, 4986113);

	mpz_t n;
	mpz_init(n);
	mpz_mul( n, P, r_star);

	mpz_t( nm1 );
	mpz_init( nm1 );
	mpz_sub_ui( nm1, n, 1 );

	mpz_t offset;
	mpz_init( offset );
	mpz_mul( offset, P, L );

	mpz_t base;
	mpz_init( base );
	mpz_set_ui( base, 2 );

	mpz_t result;
	mpz_init( result );

	for( int64_t i = 0; i < 100'000'000; i++ )
	{
		mpz_powm( result,  base,  nm1, n);
		if( mpz_cmp_si( result1, 1 ) == 0 )
		{
			std::cout << " i = " <<  i << std::endl;
		}
		mpz_add( n, n, offset);
		mpz_sub_ui( nm1, n, 1);
	}

	mpz_clear( P );
	mpz_clear( L );
	mpz_clear( r_star );
	mpz_clear( n );
	mpz_clear( nm1 );
	mpz_clear( offset );
	mpz_clear( base );
	mpz_clear( result );

	t2 = high_resolution_clock::now();
	ms_double = t2 - t1;
	std::cout << ms_double.count() << "ms\n";

}
