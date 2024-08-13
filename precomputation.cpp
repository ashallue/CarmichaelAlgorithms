#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <array>

int main()
{
  const uint64_t primes[167] = {
     3,   5,   7,  11,  13,  17,  19,  23,  29,
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


  std::vector< std::array<uint64_t, 3> > new_jobs, old_jobs, output_jobs;
  // The order of the 3-tuple {P, L, b }
  // P, the preproduct
  // L = CarmichaelLambda(P)
  // if n = PR is a CN from a 4-tuple, then the primes dividing R have to exceed b

  // The trivial preproduct
  old_jobs.push_back({1, 1, 1});

  // intended bound for the computation is 10^23
  // bounds testing is done with logarithms
  double bound = 23*log(10);

  uint64_t prime_count;
  std::cout << "How many primes to use? " ;
  std::cin >> prime_count ;

  std::cout << "The elimination rule is of the form P*L*f(p) > B where, " << std::endl;
  std::cout << " f(p) = C*p^n and p is the current prime. " << std::endl;
  std::cout << "The paper currently recommends n = 4 and C = 1" << std::endl;
  std::cout << "What power n of p do you want? " ;
  uint64_t p_exponent;
  std::cin >> p_exponent;
  std::cout << "What constant C do you want? " ;
  uint64_t C_constant;
  std::cin >> C_constant;
  std::cout << std::endl;

  bound = bound - log(C_constant);

  for( uint64_t i = 0; i < prime_count; i ++ )
  {
        while( !old_jobs.empty() )
    {
        uint64_t P = old_jobs.back()[0];
        uint64_t L = old_jobs.back()[1];
        uint64_t b = old_jobs.back()[2];
        old_jobs.pop_back();

        if( log( P ) + log( L ) + p_exponent*log( primes[i] ) > bound )
        {
          output_jobs.push_back( { P, L, b }  );
        }
        else // so current_preproduct is small enough to create more jobs
        {
          new_jobs.push_back( { P, L, primes[i] } );
          if( std::gcd( P, primes[i] - 1 ) == 1 )
          {
            new_jobs.push_back( { P*primes[i], L*((primes[i] - 1)/ std::gcd(L, primes[i] -1 ) ), primes[i] });
          }

        }
      }
      std::cout << " The prime " << primes[i] << " has " << new_jobs.size() << " working jobs and " << output_jobs.size() << " are output jobs" << std::endl;
      old_jobs = new_jobs;
      new_jobs.clear();
  }


  // output the two jobs lists here:
  /*
  std::ofstream output_file("output_jobs.txt");
  for( int i = 0; i < output_jobs.size(); i++ )
  {
    for( int j = 0; j < 3; j++)
    {
       output_file << output_jobs[i][j] << " ";
    }
      output_file << std::endl;
  }
  output_file.close();

  std::ofstream working_file("working_jobs.txt");
  for( int i = 0; i < old_jobs.size(); i++ )
  {
    for( int j = 0; j < 3; j++)
    {
       working_file << old_jobs[i][j] << " ";
    }
      working_file << std::endl;
  }
  working_file.close();
 */

}
