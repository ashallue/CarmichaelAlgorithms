#include "Preproduct.h"
#include "rollsieve.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include "mpi.h"


int main(int argc, char * argv[])
{
    int my_rank;            // my CPU number for this process
    int proc;               // number of CPUs that we have
    
    MPI_Init(&argc, &argv);                     // Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    // Find my rank
    MPI_Comm_size(MPI_COMM_WORLD, &proc);       // Find out the number of processes!
    
    std::string file_name = "very_small_" + std::to_string(my_rank) + ".txt";
    
    uint16_t work_count = 0;
    
    // preproduct info:
    std::vector< uint64_t > P =   {        3,        5,         7,       11,      13,        15,       17,       19,       23,       29,       31,       33,       35,       37,       41,       43,       47,       51,       53,       59 };
    std::vector< uint64_t > L =   {        2,        4,         6,       10,      12,         4,       16,       18,       22,       28,       30,       10,       12,       36,       40,       42,       46,       16,       52,       58 };
    std::vector< uint64_t > AB =  { 41666666, 25000000,  17857142, 11363636,  9615384,  8333333,  7352941,  6578947,  5434782,  4310344,  4032258,  3787878,  3571428,  3378378,  3048780,  2906976,  2659574,  2450980,  2358490,  2118644 };
    std::vector< uint64_t > ub =  { 69336128, 58480355,  52275796, 44964432, 42529038, 40548014, 38891112, 37475618, 35163390, 32548728, 31833138, 31176596, 30571072, 30010008, 29000494, 28543718, 27709842, 26965592, 26622044, 25687155 };
    
    for( int i = 0; i < 20; i++)
    {
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
                    Pq.appending( SMALL, q );
                    Pq.CN_multiples_of_P( file_name );
                }
            }
            q = r.nextprime();
        }
        
    }
    
    MPI_Finalize();

    return 0;
}
