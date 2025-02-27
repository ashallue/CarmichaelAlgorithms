#include "rollsieve.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <array>

// if the bound or elimination rule is changed drastically
// the use of uint64_t to store preproducts will fail

int main()
{
 
    uint64_t start = 3;
    Rollsieve incremental_sieve( start );
    
    uint32_t p = incremental_sieve.nextprime();
    // p has the value 3 right now
    
    std::vector< std::array<uint64_t, 3> > new_jobs, old_jobs, output_jobs;
      
    old_jobs.push_back({1, 1, 1});


    // intended bound for the computation is 10^24
    // bounds testing is done with logarithms
    // double bound = 24*log(10);

    // bound for a test run
    double bound = 24 * log(10);
    
    while( !old_jobs.empty() )
    {
        while( !old_jobs.empty() )
        {
            uint64_t P = old_jobs.back()[0];
            uint64_t L = old_jobs.back()[1];
            uint64_t b = old_jobs.back()[2];
            old_jobs.pop_back();

            if( ( log( P ) + 3*log( L ) + (4)*log( p ) > bound ) || ( log( P ) + 7*log( p ) > bound ) )
            {
              output_jobs.push_back( { P, L, b }  );
            }
            else // so current_preproduct is small enough to create more jobs
            {
                new_jobs.push_back( { P, L, p } );
              
                // admissibility check to create new preproduct
                if( std::gcd( P, p-1 ) == 1 )
                {
                    new_jobs.push_back( { P*p, L*( (p-1) / std::gcd( L, p-1 ) ), p });
                }
            }
        }
        p = incremental_sieve.nextprime();
        old_jobs = new_jobs;
        new_jobs.clear();
    }



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


}
