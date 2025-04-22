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
#include <numeric>

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
        output.flush();                    
        std::cout << "line read: " << P << " " << L << " " << prime_lower << "\n";

        // ignore the (1, 1, p) job, we'll deal with it in a separate loop
        if(P != 1){
            
            // create preproduct object
            Preproduct preprod = Preproduct();
            preprod.initializing( P, L, prime_lower );
            //std::cout << "Preprod initialized";
            //gmp_printf(" P = %Zd, L = %Zd\n", preprod.P, preprod.L);
       
            // search for all carmichael numbers that complete that preproduct
            preprod.CN_multiples_of_P( cars_file );
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
                preprod.CN_multiples_of_P( cars_file );
                
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
    std::cout << ms_double.count() << std::endl;

    //output.close();
    jobs_input.close();
}

// main for testing
int main(int argc, char* argv[]){
    std::cout << "num args = " << argc << "\n";
    std::cout << "0th arg is " << argv[0] << "\n";

    std::string filename;
    int job_num;
    if(argc >= 3){
        filename = argv[1];
        job_num = atoi(argv[2]);
        std::cout << " " << filename << "and job_num = " << job_num << "\n";
    }else{
        filename = "default.txt";
        job_num = 0;
    }
    
    // this is for the 10^18 job
    uint64_t B = 1'000'000'000'000'000'000;
    std::string outfile = "cars" + std::to_string(job_num);
    tabulate_test(B, filename, outfile);

    /*
    uint64_t P = 999919;
    uint64_t L = 55440;
    uint64_t prime_lower = 1009;
    Preproduct preprod;
    preprod.initializing(P, L, prime_lower);
    std::cout << "done initializing\n";
    */
    /*
    
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
   
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    
    /*
    
    uint64_t work_count = 0;
    const uint32_t my_rank = 3;
    const uint32_t proc = 25'000;
       
    
    std::vector< uint64_t > P =  { 2065, 2069, 2071, 2077, 2081, 2083, 2087, 2089, 2091, 2093, 2095, 2099, 2101,
        2103, 2111, 2113, 2117, 2119, 2123, 2129   } ;
    std::vector< uint64_t > L =  { 348, 2068, 108, 330, 2080, 2082, 2086, 2088, 80, 132, 836, 2098, 190, 700,
        2110, 2112, 504, 324, 960, 2128   };
    std::vector< uint64_t > AB = {  60532, 60415, 60357, 60182, 60067, 60009, 59894, 59837, 59780, 59722, 59665,
        59552, 59495, 59438, 59213, 59157, 59045, 58990, 58878, 58713  } ;
    std::vector< uint64_t > ub = { 7852840, 7847776, 7845249, 7837687, 7832662, 7830154, 7825149, 7822651,
        7820156, 7817664, 7815176, 7810208, 7807729, 7805253, 7795381, 7792921, 7788009,
        7785558, 7780666, 7773350 };
    
     
        std::vector< uint64_t > P =  { 1997, 1999, 2001, 2003, 2011, 2017, 2021, 2027, 2029, 2031, 2033, 2039, 2045,
            2047, 2049, 2051, 2053, 2055, 2059, 2063  } ;
        std::vector< uint64_t > L =  { 1996, 1998, 308, 2002, 2010, 2016, 966, 2026, 2028, 676, 954, 2038, 408, 88,
            682, 876, 2052, 136, 140, 2062  };
        std::vector< uint64_t > AB = {  62593, 62531, 62468, 62406, 62158, 61973, 61850, 61667, 61606, 61546, 61485,
            61304, 61124, 61064, 61005, 60945, 60886, 60827, 60709, 60591  } ;
        std::vector< uint64_t > ub = { 7940979, 7938330, 7935684, 7933042, 7922509, 7914645, 7909420, 7901608,
            7899011, 7896418, 7893827, 7886077, 7878357, 7875790, 7873227, 7870667, 7868110,
            7865557, 7860460, 7855377 };
        
        std::vector< uint64_t > P =  { 1941, 1943, 1945, 1947, 1949, 1951, 1955, 1957, 1959, 1961, 1963, 1969, 1973,
            1977, 1979, 1981, 1985, 1987, 1991, 1993  } ;
        std::vector< uint64_t > L =  { 646, 924, 388, 290, 1948, 1950, 176, 306, 652, 468, 300, 890, 1972, 658, 1978,
            282, 396, 1986, 180, 1992  };
        std::vector< uint64_t > AB = {  64399, 64333, 64267, 64201, 64135, 64069, 63938, 63873, 63808, 63742, 63678,
            63484, 63355, 63227, 63163, 63099, 62972, 62908, 62782, 62719  } ;
        std::vector< uint64_t > ub = { 8016625, 8013873, 8011126, 8008382, 8005641, 8002905, 7997443, 7994718,
            7991996, 7989278, 7986564, 7978444, 7973048, 7967667, 7964982, 7962301, 7956949,
            7954279, 7948948, 7946288 };
        
        std::vector< uint64_t > P =  { 1883, 1885, 1889, 1891, 1895, 1897, 1901, 1903, 1907, 1909, 1913, 1915, 1919,
            1921, 1923, 1927, 1931, 1933, 1937, 1939  } ;
        std::vector< uint64_t > L =  { 804, 84, 1888, 60, 756, 270, 1900, 860, 1906, 902, 1912, 764, 900, 112, 640,
            920, 1930, 1932, 444, 276 };
        std::vector< uint64_t > AB = {  66383, 66312, 66172, 66102, 65963, 65893, 65754, 65685, 65547, 65479, 65342,
            65274, 65138, 65070, 65002, 64867, 64733, 64666, 64532, 64466  } ;
        std::vector< uint64_t > ub = { 8098103, 8095238, 8089520, 8086667, 8080973, 8078132, 8072463, 8069634,
            8063988, 8061170, 8055548, 8052743, 8047144, 8044350, 8041560, 8035992, 8030440,
            8027669, 8022139, 8019380 };
        
        
        std::vector< uint64_t > P =  { 1817, 1819, 1823, 1829, 1831, 1835, 1837, 1841, 1843, 1847, 1851, 1853, 1855,
            1861, 1865, 1867, 1871, 1873, 1877, 1879 } ;
        std::vector< uint64_t > L =  { 858, 848, 1822, 870, 1830, 732, 830, 786, 288, 1846, 616, 432, 156, 1860, 372,
            1866, 1870, 1872, 1876, 1878 };
        std::vector< uint64_t > AB = {  68794, 68719, 68568, 68343, 68268, 68119, 68045, 67897, 67824, 67677, 67531,
            67458, 67385, 67168, 67024, 66952, 66809, 66737, 66595, 66524 } ;
        std::vector< uint64_t > ub = { 8194990, 8191986, 8185990, 8177028, 8174050, 8168106, 8165141, 8159223,
            8156271, 8150379, 8144503, 8141572, 8138645, 8129889, 8124073, 8121171, 8115379,
            8112490, 8106723, 8103846 };
        
        std::vector< uint64_t > P =  { 1753, 1757, 1759, 1761, 1763, 1765, 1769, 1777, 1779, 1781, 1783, 1787, 1789,
            1793, 1795, 1797, 1799, 1801, 1807, 1811 } ;
        std::vector< uint64_t > L =  { 1752, 750, 1758, 586, 840, 352, 420, 1776, 592, 408, 1782, 1786, 1788, 810,
            716, 598, 768, 1800, 276, 1810 };
        std::vector< uint64_t > AB = {71306, 71143, 71063, 70982, 70901, 70821, 70661, 70343, 70264, 70185, 70106,
            69949, 69871, 69715, 69637, 69560, 69483, 69405, 69175, 69022  } ;
        std::vector< uint64_t > ub = { 8293530, 8287232, 8284090, 8280952, 8277820, 8274692, 8268451, 8256024,
            8252929, 8249838, 8246753, 8240595, 8237523, 8231393, 8228334, 8225281, 8222231,
            8219187, 8210079, 8204030  };
   
 
    
    
    for( int i = 0; i < 20; i++)
    {
        t1 = high_resolution_clock::now();
       
        Rollsieve r( AB[i] );
        uint64_t q = r.nextprime();
        Preproduct small_P = Preproduct();
        small_P.initializing( P[i], L[i], AB[i] );
        Preproduct Pq = Preproduct();
       
        while( q < ub[i] )
        {
            if( small_P.is_admissible_modchecks( q ) )
            {
                work_count++;
                work_count = work_count % proc;
                if( work_count == my_rank )
                {
                    Pq.appending( small_P, q );
                    Pq.CN_multiples_of_P( "new_CNs_block_timing.txt" );
                }
            }
            q = r.nextprime();
        }
        
        t2 = high_resolution_clock::now();
        ms_double = t2 - t1;
        std::cout << P[i] << std::endl;
        std::cout << ms_double.count() << std::endl;
   }
    */
 
    /*
    uint64_t P = 921041 ;
    uint64_t L = 360 ;
    uint64_t AB = 73 ;
    
    t1 = high_resolution_clock::now();
   
    Rollsieve r( AB );
    uint64_t q = r.nextprime();
    Preproduct small_P = Preproduct();
    small_P.initializing( P, L, AB );
    Preproduct Pq = Preproduct();
   
    while( q < 10278 )
    {
        if( small_P.is_admissible_modchecks( q ) )
        {
            std::cout << " being called with q = " << q << std::endl;  //stalled on q = 251
            Pq.appending( small_P, q );
            Pq.CN_multiples_of_P( "new_CNs_block_timing.txt" );
        }
        q = r.nextprime();
    }
    
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << P << std::endl;
    std::cout << ms_double.count() << std::endl;
     
     uint64_t P = 231181291 ;
     uint64_t L = 9000 ;
     uint64_t AB = 251 ;
     
     t1 = high_resolution_clock::now();
    
     Rollsieve r( AB );
     uint64_t q = r.nextprime();
     Preproduct small_P = Preproduct();
     small_P.initializing( P, L, AB );
     Preproduct Pq = Preproduct();
    
     while( q < 65770 )
     {
         if( small_P.is_admissible_modchecks( q ) )
         {
             std::cout << " being called with q = " << q << std::endl;  //stalled on q = 811
             Pq.appending( small_P, q );
             Pq.CN_multiples_of_P( "new_CNs_block_timing.txt" );
         }
         q = r.nextprime();
     }
     
     t2 = high_resolution_clock::now();
     ms_double = t2 - t1;
     std::cout << P << std::endl;
     std::cout << ms_double.count() << std::endl;
     
     
     */

    /*
    uint64_t P = 187488027001 ;
    uint64_t L = 81000 ;
    uint64_t AB = 811 ;
    
    t1 = high_resolution_clock::now();
   
    Preproduct small_P = Preproduct();
    small_P.initializing( P, L, AB );
   
    std::cout << "calling on suspected problem case: " << std::endl;
    small_P.CN_search( "new_CNs_block_timing.txt" );

    
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << P << std::endl;
    std::cout << ms_double.count() << std::endl;
    */
    
    
    /*
    std::cout << "starting timing test for job (" << P << ", " << L << ", " << AB << ")\n";
    job_timing(P, L, AB, "single_job.txt");
    


    Preproduct PP = Preproduct();
    PP.initializing( 41, 40, 100'000/41 );

    PP.CN_multiples_of_P( "CN_divis_7.txt" );
    PP.initializing( 43, 42, 100'000/43 );

    PP.CN_multiples_of_P( "CN_divis_7.txt" );
    PP.initializing( 47, 46, 100'000/47 );

    PP.CN_multiples_of_P( "CN_divis_7.txt" );
    PP.initializing( 51, 16, 100'000/51 );

    PP.CN_multiples_of_P( "CN_divis_7.txt" );
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
