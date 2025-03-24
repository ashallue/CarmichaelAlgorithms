/*  Test functions for tabulating carmichaels project.
 *  Andrew Shallue and Jonathan Webster, January 2025

First several carmichaels with 3 prime factors
561 3 11 17
1105 5 13 17
1729 7 13 19
2465 5 17 29
2821 7 13 31
6601 7 23 41
8911 7 19 67
10585 5 29 73
15841 7 31 73
29341 13 37 61
with 4 prime factors:
41041 7 11 13 41
62745 3 5 47 89
63973 7 13 19 37
75361 11 13 17 31
101101 7 11 13 101
126217 7 13 19 73
172081 7 13 31 61
188461 7 13 19 109
278545 5 17 29 113
340561 13 17 23 67

Carmichaels with 3 prime factors whose first factor is 3 or 5 or 7
561 3 11 17
ashallue@hyperion:~/tabulate_car/datafiles_22$ awk 'NF == 4 && $2 == 5' cars_table_10to22.txt 
1105 5 13 17
2465 5 17 29
10585 5 29 73
ashallue@hyperion:~/tabulate_car/datafiles_22$ awk 'NF == 4 && $2 == 7' cars_table_10to22.txt 
1729 7 13 19
2821 7 13 31
6601 7 23 41
8911 7 19 67
15841 7 31 73
52633 7 73 103
 */

#include "test_functions.h"
#include <stdio.h>

// call CN_factorization on all squarefree multiples of 7 with 4 prime factors
bool test_factor(){
    // file for printing carmichael numbers
    std::string cars_file = "cars.txt";
    
    // primes for building numbers to factor
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 
    uint64_t num_small_primes = 20;
    uint64_t* small_primes = new uint64_t[ num_small_primes ];
    small_primes[0] = 5;  small_primes[1] = 7; small_primes[2] = 11; small_primes[3] = 13; 
    small_primes[4] = 17;  small_primes[5] = 19; small_primes[6] = 23; small_primes[7] = 29; 
    small_primes[8] = 31;  small_primes[9] = 37; small_primes[10] = 41; small_primes[11] = 43; 
    small_primes[12] = 47;  small_primes[13] = 53; small_primes[14] = 59; small_primes[15] = 61;
    small_primes[16] = 67;  small_primes[17] = 71; small_primes[18] = 73; small_primes[19] = 79; 

    // initialize preproduct with value 7.
    Preproduct PPP;
    uint64_t p = 7;
    uint64_t pl = 6;
    uint64_t bound = 1000;
    PPP.initializing(p, pl, bound);

    // initialize n, r, and the vector that will store prime factors of r.  Note n = P * r
    mpz_t n;
    mpz_init(n);
    mpz_t r;
    mpz_init(r);
    std::vector<uint64_t> r_primes;

    // bound on the index of the prime factors
    uint64_t factor_bound = 19;

    // tracking variables
    uint64_t count_factored = 0;
    uint64_t count_not_carmichael = 0;
    bool factor_fail = true;
    
    // we will factor multiples of 3 where we know the factors
    for(uint64_t i1 = 0; i1 < factor_bound; ++i1){
        for(uint64_t i2 = i1 + 1; i2 < factor_bound; ++i2){
            for(uint64_t i3 = i2 + 1; i3 < factor_bound; ++i3){

                // test for 41041 7 11 13 41
                if( small_primes[i1] == 11 && small_primes[i2] == 13 && small_primes[i3] == 41 ){
                    std::cout << "41041 = 7 * 11 * 13 * 41 found\n";
                }
                // reset r, n
                mpz_set_ui(r, 1);
                mpz_set_ui(n, p);
                
                // form r and n
                mpz_mul_ui( r, r, small_primes[i1] );
                mpz_mul_ui( r, r, small_primes[i2] );
                mpz_mul_ui( r, r, small_primes[i3] );
                mpz_mul(n, n, r);

                // testing
                //gmp_printf ("n = %Zd ", n );  std::cout << " ";
                //gmp_printf ("r = %Zd ", r );  std::cout << "\n";
                
                // factor
                
                bool is_factored = PPP.CN_factorization(n, r, cars_file);

                /* testing
                if( small_primes[i1] == 11 && small_primes[i2] == 13 && small_primes[i3] == 41 ){
                    std::cout << "41041, CN_factorization returned " << is_factored << "\n";
                    std::cout << "r_primes: ";
                    for(int i = 0; i < r_primes.size(); ++i){
                        std::cout << r_primes.at(i) << " ";
                    }
                    std::cout << "\n";
                }
                */
                
                // check whether r is fully factored
                bool factors_match = r_primes.size() == 3 && r_primes.at(0) == small_primes[i1] &&
                    r_primes.at(1) == small_primes[i2] && r_primes.at(2) == small_primes[i3];

                // update tracking variables
                if(is_factored){
                    count_factored++;
                    if( !factors_match ) factor_fail = false;
                }else{
                    count_not_carmichael++;
                }
                
           
                    //gmp_printf ("%Zd = ", n );
                    //std::cout << " is carmichael \n";
            }
        }
    }
    std::cout << count_factored << " many numbers factored, " << count_not_carmichael << " failed a fermat test\n";
    
    delete[] small_primes;
    
    return factor_fail;
}

// function that tabulates all Carmichaels up to a given bound.  Only large case, no small case
// need to read in jobs from output_jobs.txt, create a preproduct object for each one, call search function, assemble the output
// Currently the method variable would be CN_search or CN_search_no_wheel
void tabulate_test(uint64_t bound, std::string jobs_file, std::string cars_file){
    //setup jobs as an input file, setup file for writing carmichael numbers as output
    std::ofstream output;
    output.open(cars_file, std::ios_base::app);
    std::ifstream jobs_input;
    jobs_input.open(jobs_file);

    //std::string input_line = jobs_input.getline();
    uint64_t P, L, prime_lower;
    uint64_t prime;
    uint64_t third_bound = pow(bound, 1.0/3);

    std::cout << "tabulating up to " << bound << " with primes up to " << third_bound << "\n";
    
    // test that I can read and process a job triple
    while(jobs_input >> P >> L >> prime_lower){
        std::cout << "line read: " << P << " " << L << " " << prime_lower << "\n";

        // ignore the (1, 1, p) job, we'll deal with it in a separate loop
        if(P != 1){
            
            // create preproduct object
            Preproduct preprod = Preproduct();
            preprod.initializing( P, L, prime_lower );

            // search for all carmichael numbers that complete that preproduct
            preprod.complete_tabulation( cars_file );
        }else{
            // time to deal with the (1, 1, p) job.  We replace with (p, p-1, p) for all 
            // primes p satisfying prime_lower < p < bound^(1/3)

            // set up a rollsieve to generate primes
            Rollsieve prime_gen = Rollsieve(prime_lower + 1);
            prime = prime_gen.nextprime();

            // loop until bound reached
            while ( prime < third_bound ){
                std::cout << "Running job for prime " << prime << "\n";

                // for (p, p-1, p) job
                P = prime;
                L = prime - 1;
                prime_lower = prime;

                // create preproduct object and run job
                Preproduct preprod = Preproduct();
                preprod.initializing( P, L, prime_lower );
                preprod.complete_tabulation( cars_file );
                
                // next prime
                prime = prime_gen.nextprime();
            }
        }
    } 
    //output.close();
    jobs_input.close();
}

// function that prints wall time for one particular job
void job_timing(uint64_t P, uint64_t L, uint64_t prime_lower, std::string cars_file){
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    //setup jobs as an input file, setup file for writing carmichael numbers as output
    std::ofstream output;
    output.open(cars_file, std::ios_base::app);
    std::ifstream jobs_input;
    jobs_input.open(cars_file);
    
    auto t1 = high_resolution_clock::now();

    // create preproduct object and run job
    Preproduct preprod = Preproduct();
    preprod.initializing( P, L, prime_lower );

    preprod.CN_multiples_of_P( cars_file );
    
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout << ms_double.count() << "ms\n";

    //output.close();
    jobs_input.close();
}

// main for testing
int main(){
    
   
    uint64_t work_count = 0;
    const uint16_t my_proc_number = 3;
    const uint16_t total_procs = 25'000;
    
    uint64_t P;
    uint64_t L;
    uint64_t AB;
    
    Rollsieve r(41'666'667);
    uint64_t q = r.nextprime();

    while( q < 69'336'127 )
    {
        if( 1 != q % 3 )
        {
            work_count++;
            if( my_proc_number == work_count % total_procs )
            {
                P = 3*q;
                L = q-1;
                AB = q;
                std::cout << "starting timing test for job (" << P << ", " << L << ", " << AB << ")\n";
                job_timing(P, L, AB, "single_job.txt");
            }
        }
        q = r.nextprime();
    }
    
     
    /*
    uint64_t P = 3 ;
    uint64_t L = 2 ;
    uint64_t AB = 41'666'667;
    
    std::cout << "starting timing test for job (" << P << ", " << L << ", " << AB << ")\n";
    job_timing(P, L, AB, "single_job.txt");
    */
    
     /*
    std::cout << "\nTabulating up to 10^15\n";
    uint64_t upper = 1000000000000000;
    tabulate_test(upper, "output_jobs.txt", "small_tabulation.txt");
    */
    
      
      
    //Preproduct preprod = Preproduct();
    //preprod.initializing( 1, 1, 97 );

    // search for all carmichael numbers that complete that preproduct
    //preprod.CN_search( "cars.txt" );
    
    /*
    // testing file writing.  Source: https://cplusplus.com/reference/cstdio/fprintf/
    FILE * pFile;
   int m;
   char name [100];

   pFile = fopen ("myfile.txt","w");
   for (m=0 ; m<3 ; m++)
   {
     puts ("please, enter a name: ");
     std::cin >> name;
     fprintf (pFile, "Name %d [%-10.10s]\n",m+1,name);
     gmp_fprintf(pFile, "%Zd = ", n);
   }
   fclose (pFile);
    */
    




}
