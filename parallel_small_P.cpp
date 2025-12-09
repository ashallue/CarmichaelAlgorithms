#include "Preproduct.h"
#include "rollsieve.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include "mpi.h"

// compiled and run with:
// mpic++ -o parallel parallel_small_P.cpp Preproduct.o rollsieve.o -lgmp -O3
// mpirun --hostfile hostfile.txt -n 156 ./parallel &

int main(int argc, char * argv[])
{
    int my_rank;            // my CPU number for this process
    int proc;               // number of CPUs that we have

    MPI_Init(&argc, &argv);                     // Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    // Find my rank
    MPI_Comm_size(MPI_COMM_WORLD, &proc);       // Find out the number of processes!

    std::string file_name = "nbd4_53_59_very_small_" + std::to_string(my_rank) + ".txt";
    
    uint16_t work_count = 0;

    std::vector< uint64_t > P =  { 53, 59 };
    std::vector< uint64_t > L =  { 52, 58 };
    std::vector< uint64_t > AB =  { 2358490, 2118644 };
    std::vector< uint64_t > ub =  { 26622044, 25687156 };


    uint16_t num_cases = P.size();
    
    for( int i = 0; i < num_cases; i++ )
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
                    Pq.appending( small_P, q );
                    Pq.CN_multiples_of_P( file_name );
                }
            }
            q = r.nextprime();
        }
    }

    MPI_Finalize();

    return 0;
}

/*
These have been done:

(0) old bigdawg     1 < P < 60          started 3-25: ended April 7th
new bigdawg     60 < P < 114        started 3-26: less than 48 hours
xeon phi        P = 115             started 3-27: 20 hours
xeon phi        P = 119             started 3-28: 22 hours
new bigdawg     120 < P < 178       started 8:21am 3-29: approx 24 hours
xeon phi        178 < P < 186       started 8:24am 3-29: approx 61 hours
nbd node-4      P = 187             started 2:27pm 3-30: approx 8 hours
xeon phi        P = 191             started 8:08am 4-1: 21 hours
nbd             192 < P < 240       started 7:40am 3-31: approx 25 hours
nbd             240 < P < 304       started 2:00pm 4-1: approx 27 hours
xeon phi        P = 307             started 9:25am 4-2: 17 hours
xeon phi        310 < P < 320       started 7:45am 4-3: 62 hours
nbd             320 < P < 368       started 7:45am 4-3: about 16 hours
nbd node-0      370 < P < 430       started 6:00am 4-4: less than 4 days
nbd node-1      430 < P < 486       started 6:00am 4-4: a little over 3 days
nbd node-2      486 < P < 548       started 6:00am 4-4:  under 3 days
nbd node-3      548 < P < 602       started 6:00am 4-4:  under 3 days
nbd node-3      602 < P < 674       started 10:40am 4-7: a little over 2 days
nbd node-2      674 < P < 732       started 10:40am 4-7: a little over 2 days
xeon phi        732 < P < 794       started 10:40am 4-7: 6 days, 15 hours
nbd node-1      794 < P < 858       started 3:10pm 4-7: a little under 2 days
(1) obd             858 < P < 920       started 7:00am 4-8: less than 60 hours
nbd node-0      920 < P < 990       started 7:00am 4-8: a little over 2 days
xeon phi        P = 991             started 9:25am 4-14: 8.4 hours
nbd node-0      P = 995             started 10:35am 4-14: under 3 hours
nbd node-1      P = 997             started 11:10am 4-14: under 3 hours
nbd node-2      P = 1001            started 12:15am 4-14: under 3 hours
nbd node-3      P = 1003            started 12:15am 4-14: under 3 hours
nbd node-4      P = 1007            started 12:15am 4-14: under 3 hours
nbd node-0      P = 1009            started 1:00pm, 4-14: under 3 hours
nbd node-1      P = 1013            started 1:10pm, 4-14: under 3 hours
nbd node-3      P = 1019            started 2:20pm, 4-14: under 3 hours
xeon phi        P = 1021            started 9:10am, 4-15: 8 hours
(2) obd             1022 < P < 1060     started 11:30am 4-15: less than 48 hours
nbd node-0      1060 < P < 1124     started 7:00am 4-10: about 41 hours
(3) obd             1124 < P < 1184     started 6:45am 4-17: 60 hours
nbd node-3      1184 < P < 1246     started 7:00am 4-10: about 38 hours
(4) obd             1246 < P < 1362     started 2:00pm 4-19: about 84 hours ~ 2.1 hours per P
nbd node-0      1362 < P < 1430     started 8:00am 4-12: less than 36 hours
nbd node-3      1430 < P < 1494     started 8:00am 4-12: less than 36 hours
xeon phi        P = 1495            started 8:00pm 4-16: about 9 hours
xeon phi        P = 1499            started 2:00pm 4-19: about 5 hours
xeon phi        P = 1501            started 7:00am 4-21: less than 8 hours
xeon phi        1506 < P < 1528     started 3:00pm 4-21: 29 hours ~ 3.6 hours per P
xeon phi        1528 < P < 1554     started 7:40am 4-23: 50 hours ~ 5 hours per P
nbd node-0      1554 < P < 1750     started 9:20am 4-23: 4-26 9:00pm
nbd node-1      1750 < P < 1940     started 9:20am 4-23: 4-26 4:00pm
nbd node-2      1940 < P < 2130     started 9:20am 4-23: 72 hours
nbd node-3      2130 < P < 2334     started 9:20am 4-23: 72 hours
nbd node-4      2334 < P < 2514     started 9:20am 4-23: 65 hours ~  1.1 hours per P
(5) obd             2514 < P < 2712     started 9:00am 4-24: 86 hours ~ 1.4 hours per P
xeon phi        2712 < P < 2768     started 10:25am 4-25: 62 hours ~ 3 hours per P
nbd node-4      2768 < P < 3074     started 6:30am 4-26: 4-30 3:00am
nbd node-2      3074 < P < 3406     started 3:45pm 4-27: about 88 hours
nbd node-3      3406 < P < 3728     started 1:00pm 4-26: less than 80 hours ~ .8 hours per P
nbd node-0      3728 < P < 4062     started 3:45pm 4-27: 71 hours
nbd node-1      4062 < P < 4380     started 3:45pm 4-27: 70 hours ~ .7 hours per P
(6) obd             4380 < P < 4704     started 4:00am 4-28: 99 hours ~ 1 hour per P
xeon phi        4704 < P < 4766     started 7:00am 4-28: 35 hours ~ 1.75 hours per P
nbd node-3      4766 < P < 5100     started 8:00am 4-30: 60 hours
xeon phi        5100 < P < 5168     started 8:00pm 4-29: 35 hours
nbd node-4      5168 < P < 5816     started 8:00am 4-30: 113 hours
nbd node-1      5816 < P < 6458     started 1:30pm 4-30: 102 hours
nbd node-0      6458 < P < 7134     started 3:00pm 4-30: 96 hours ~ .48 hours per P
xeon phi        7134 < P < 7250     started 8:00am  5-1: 48 hours ~ 1.2 hours per P
nbd node-2      7250 < P < 7896     started 8:00am  5-1: ended 4-4 10:00pm
(7) obd             7896 < P < 8532     started 10:00am 5-2: 121 hours ~ 1.65 P per hour
nbd node-3      8532 < P < 9170     started 10:00am 5-2: finished
xeon phi        9170 < P < 9504     started 9:00am 5-3: 100 hours  ~ 1 P per hour
nbd node-0      9504 < P < 10308    started 3:00pm 5-4: 89 hours
nbd node-1      10308 < P < 11162   started 7:00am 5-5: 81 hours ~ 3.1 P per hour
nbd node-2      11162 < P < 11988   started 7:00am 5-5: 74 hours ~ 3.4 P per hour
nbd node-4      11988 < P < 12804   started 7:00am 5-5: 70 hours ~ 3.6 P per hour
nbd node-3      12804 < P < 13592   started 7:00am 5-6: 70.5 hours
(8) obd             13592 < P < 14894   started 11:40am 5-7:  finished 5-13
xeon phi        14894 < P < 15572   started 2:00pm 5-7: 144 hours
nbd node-4      15572 < P < 17224   started 6:30am 5-8: 5-12 at 10 pm
nbd node-0      17224 < P < 18882   started 8:00am 5-8: 97 hours ~ 5.2 P per hour
nbd node-2      18882 < P < 20522   started 9:20am 5-8: 93 hours ~ 5.4 P per hour
nbd node-1      20522 < P < 22164   started 4:10pm 5-8: 84 hours ~ 6.0 P per hour
nbd node-3      22164 < P < 23802   started 7:50am 5-9: 79 hours ~ 6.3 P per hour
nbd node-0      23802 < P < 25458   started 10:10am 5-19: 75 hours ~ 6.7 P per hour
nbd node-1      25458 < P < 27484   started 10:10am 5-19: 81 hours
nbd node-2      27484 < P < 29788   started 10:10am 5-19: 90 hours
nbd node-3      29788 < P < 32454   started 10:10am 5-19: 94 hours
nbd node-4      32454 < P < 35442   started 10:10am 5-19: 98 hours
xeon phi        35442 < P < 36090   started 10:10am 5-19: 64.5 hours ~ 3.1 P per hour
xeon phi        36090 < P < 37748   started 6:50am 5-22: 141 hours ~ 3.5 P per hour
nbd node-0      37748 < P < 41090   started 1:40pm 5-22: 98 hours
nbd node-1      41090 < P < 44354   started 8:00am 5-23: 91 hours
nbd node-2      44354 < P < 47688   started 8:00am 5-23: 85 hours
nbd node-3      47688 < P < 51042   started 8:40am 5-23: 78 hours ~ 12.8 P per hour
nbd node-4      51042 < P < 54386   started 2:00pm 5-23: 74 hours ~ 13.5 P per hour
nbd node-3      54386 < P < 56060   started 3:30pm 5-26: 36 hours ~ 13.9 P per hour
nbd node-0      56060 < P < 57728   started 3:50pm 5-26: 33 hours ~ 15.0 P per hour
nbd node-4      57728 < P < 59414   started 4:10pm 5-26: 33 hours ~ 15.0 P per hour
nbd node-1      59414 < P < 60410   started 6:30am 5-27: 21 hours ~
(9) obd             60410 < P < 67048   started 6:30pm 5-13: 184 hours ~ 10.8 P per hour
(10) obd             67048 < P < 73722   started 11:10am 5-21: 167 hours ~ 12 P per hour
nbd node-2      73722 < P < 74752   started 11:20am 5-27: 15 hours ~ 20 P per hour
nbd node-0      74752 < P < 76428   started 7:40am 5-28: 26 hours ~ 19 P per hour
nbd node-1      76428 < P < 78132   started 7:40am 5-28: 26 hours
nbd node-2      78132 < P < 79784   started 7:45am 5-28: 26 hours
nbd node-3      79784 < P < 81464   started 7:45am 5-28: 24 hours
nbd node-4      81464 < P < 83138   started 7:50am 5-28: 24.5 hours
(11) obd             83138 < P < 89840   started 10:55am 5-28: 141 hours  ~ 14 P per hour
xeon phi        89840 < P < 91549   started 8:30am 5-29: 64.5 hours
 
Old bigdawg used about 220k core hours:  XEON E5-2630 v2
We think this is about 110k core hours on EPYC 9734. 
Repair computation by regime number from above:

 (0)
 1 < P < 30: 2 preproducts per node.  Started at 7:00am 8-14.
 30< P < 60:
 
 (1)
 858 < P < 920:  4 preproduct per node.  Started at 12:40pm 8-13.  finished 10 hours later
 
 (2)
 1022 < P < 1060: 2 preproducts per node.   Started at 7:00am 8-13.  finished within 5 hours.
 
 (3)
 1124 < P < 1184:  4 preproducts per node.  Started at 8:30am 8-12  finished 9-10 hours later
 
 (4)
 nbd node-2      [ 1247, 1277 ]    started 3:10pm 8-10   23 hours
 nbd node-3      [ 1279, 1303 ]    started 3:10pm 8-10:  15 hours
 nbd node-0      [ 1307, 1335 ]    started 7:40am 8-11:  14 hours
 nbd node-3      [ 1337, 1361 ]    started 7:40am 8-11:  15 hours
 
 (5)
 nbd node-1      [ 2515, 2711 ]    started 9:00am 8-9:  7:00pm 8-11
 
 (6)
 nbd node-0      [ 4381, 4703 ]    started 10:24am 8-8:  5:00am 8-11
 
 (7)
 nbd node-2      [ 7897, 8221 ]    started 8:20am 8-8:   3:00am 8-10
 nbd node-3      [ 8223, 8531 ]    started 8:20am 8-8:   11:00pm 8-9
 
 (8)
 nbd node-3      [ 13593, 13865 ]  started 7:00am 8-7:   21 hours
 nbd ndoe-2      [ 13871, 14131 ]  started 10:30am 8-7:  21 hours
 nbd node-4      [ 14137, 14893 ]  started 8:15am 8-8:  5:30pm 8-10
 
 (9)
 nbd node-2      [ 60413, 61085 ]  started 2:00pm 8-4:  13 hours
 nbd node-3      [ 61087, 61763 ]  started 2:10pm 8-4:  13 hours
 nbd node-2      [ 61769, 62085 ]  started 7:24am 8-5:  6.75 hours
 nbd node-3      [ 62087, 62401 ]  started 1:25pm 8-5:  7 hours
 nbd node-2      [ 62405, 63393 ]  started 2:15pm 8-5:  19 hours
 nbd node-3      [ 63397, 64081 ]  started 8:00am 8-6:  12 hours
 nbd node-2      [ 64083, 65399 ]  started 9:10am 8-6:  25 hours
 nbd node-4      [ 65401, 66737 ]  started 7:40pm 8-6:  24 hours
 nbd node-3      [ 66739, 67047 ]  started 8:00pm 8-6:  6 hours
 
 (10)
 nbd node-1      [ 67048, 73722 ]  started 11:30am 8-4:  117 hours

 (11)
 nbd node-0      [ 83138, 89840 ]  started 9:50am 8-4:  95 hours (2 slow threads)

 
magma script generated the inputs:

 num_preprods:= 400; P:=[]; L:=[]; AB:= []; ub:= []; count:= 0; i:= 13593;
 while count lt num_preprods do
     if GCD( i, EulerPhi( i ) ) eq 1 then
         i_f:= Factorization( i ) ;
         if i*i_f[ #i_f ][1]^3 lt 10^24 then
             count:= count + 1;
             P[count]:= i;
             L[count]:= CarmichaelLambda(i);
             AB[count]:= Max( Floor(125000000/i), i_f[ #i_f ][1] ) ;
             ub[count]:= Ceiling( (10^24/i)^(1/3) + 1);
         end if;
     end if;
     i:= i + 2;
 end while;
 printf "std::vector< uint64_t > P =  { %o", P[1] ;
 for i := 2 to num_preprods do printf ", %o", P[i]; end for;
 printf " };  \n" ;
 printf "std::vector< uint64_t > L =  { %o", L[1] ;
 for i := 2 to num_preprods do printf ", %o", L[i]; end for;
 printf " };  \n" ;
 printf "std::vector< uint64_t > AB =  { %o", AB[1] ;
 for i := 2 to num_preprods do printf ", %o", AB[i]; end for;
 printf " };  \n" ;
 printf "std::vector< uint64_t > ub =  { %o", ub[1] ;
 for i := 2 to num_preprods do printf ", %o", ub[i]; end for;
 printf " };  \n" ;
   


 
 */
