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

    std::string file_name = "nbd_node3_551_601_small_" + std::to_string(my_rank) + ".txt";

    uint16_t work_count = 0;

     std::vector< uint64_t > P =  {      607,      611,      613,      617,      619,      623,      629,      631,      635,      641,      643,      647,      649,      653,      659,      661,      665,      667,      671,      673 };
     std::vector< uint64_t > L =  {      606,      276,      612,      616,      618,      264,      144,      630,      252,      640,      642,      646,      290,      652,      658,      660,       36,      308,       60,      672 };
     std::vector< uint64_t > AB = {   205930,   204582,   203915,   202593,   201938,   200642,   198728,   198098,   196850,   195007,   194401,   193199,   192604,   191424,   189681,   189107,   187969,   187406,   186289,   185735 };
     std::vector< uint64_t > ub = { 11810560, 11784730, 11771900, 11746406, 11733741, 11708575, 11671227, 11658883, 11634351, 11597936, 11585899, 11561973, 11550084, 11526452, 11491364, 11479763, 11456699, 11445237, 11422449, 11411123 };
    /*
     std::vector< uint64_t > P =  {      677,      679,      681,      683,      685,      691,      695,      697,      699,      701,      703,      705,      707,      709,      713,      717,      719,      721,      727,      731 };
     std::vector< uint64_t > L =  {      676,       96,      226,      682,      136,      690,      276,       80,      232,      700,       36,       92,      300,      708,      330,      238,      718,      102,      726,      336 };
     std::vector< uint64_t > AB = {   184638,   184094,   183553,   183016,   182481,   180897,   179856,   179340,   178826,   178316,   177809,   177304,   176803,   176304,   175315,   174337,   173852,   173370,   171939,   170998 };
     std::vector< uint64_t > ub = { 11388604, 11377411, 11366263, 11355157, 11344095, 11311166, 11289424, 11278616, 11267849, 11257122, 11246437, 11235792, 11225187, 11214622, 11193611, 11172757, 11162388, 11152057, 11121292, 11100970 };
     
     std::vector< uint64_t > P =  {      733,      739,      743,      745,      749,      751,      753,      757,      761,      763,      767,      769,      771,      773,      779,      781,      785,      787,      789,      793 };
     std::vector< uint64_t > L =  {      732,      738,      742,      148,      318,      750,      250,      756,      760,      108,      348,      768,      256,      772,      360,       70,      156,      786,      262,       60 };
     std::vector< uint64_t > AB = {   170532,   169147,   168236,   167785,   166889,   166444,   166002,   165125,   164257,   163826,   162972,   162548,   162127,   161707,   160462,   160051,   159235,   158831,   158428,   157629 };
     std::vector< uint64_t > ub = { 11090865, 11060767, 11040882, 11030994, 11011322, 11001538, 10991789, 10972395, 10953137, 10943558, 10924501, 10915022, 10905576, 10896162, 10868115, 10858830, 10840355, 10831164, 10822005, 10803778 };
     
     std::vector< uint64_t > P =  {      795,      797,      799,      803,      805,      807,      809,      811,      815,      817,      821,      823,      827, 829, 835,         839, 843, 851, 853, 857  } ;
     std::vector< uint64_t > L =  {       52,      796,      368,      360,      132,      268,      808,      810,      324,      126,      820,      822,      826, 828, 332, 838,         280, 396, 852, 856 };
     std::vector< uint64_t > AB = {   157232,   156838,   156445,   155666,   155279,   154894,   154511,   154130,   153374,   152998,   152253,   151883,   151148, 150784, 149700, 148986, 148279, 146886, 146541,         145857 } ;
     std::vector< uint64_t > ub = { 10794711, 10785674, 10776667, 10758743, 10749826, 10740938, 10732080, 10723250, 10705678, 10696935, 10679535, 10670877, 10653645, 10645071, 10619512, 10602609, 10585813, 10552537, 10544283, 10527852 };
     
     std::vector< uint64_t > P =  {   859, 863, 865, 869, 871, 877, 879, 881, 883, 885, 887, 893, 895, 899, 901,         907, 911, 913, 917, 919 } ;
     std::vector< uint64_t > L =  {  858, 862, 172, 390, 132, 876, 292, 880, 882, 116, 886, 414, 356, 420, 208,         906, 910, 410, 390, 918 };
     std::vector< uint64_t > AB = {  145518, 144843, 144508, 143843, 143513, 142531, 142207, 141884, 141562,         141242, 140924, 139977, 139664, 139043, 138734, 137816, 137211, 136911, 136314,         136017} ;
     std::vector< uint64_t > ub = { 10519675, 10503397, 10495296, 10479168, 10471141, 10447207, 10439277, 10431372, 10423490, 10415632, 10407798, 10384436, 10376695, 10361282, 10353610, 10330729, 10315587, 10308049, 10293039, 10285567 };
     
     std::vector< uint64_t > P =  {   923, 929, 933, 937, 941, 943, 947, 949, 951, 953, 957, 959, 965, 967, 971,         973, 977, 983, 985, 989 } ;
     std::vector< uint64_t > L =  { 420, 928, 310, 936, 940, 440, 946, 72, 316, 952, 140, 408, 192, 966, 970, 138,         976, 982, 196, 462 };
     std::vector< uint64_t > AB = {  135427, 134553, 133976, 133404, 132837, 132555, 131995, 131717, 131440,         131164, 130616, 130344, 129533, 129265, 128733, 128468, 127942, 127161, 126903,         126390   } ;
     std::vector< uint64_t > ub = { 10270687, 10248528, 10233861, 10219278, 10204777, 10197558, 10183180, 10176021, 10168882, 10161764, 10147586, 10140527, 10119467, 10112485, 10098580, 10091656, 10077865, 10057319, 10050507, 10036939 };
     
     std::vector< uint64_t > P =  {   991, 995, 997, 1001, 1003, 1007, 1009, 1013, 1019, 1021, 1031, 1033, 1037,         1039, 1041, 1043, 1049, 1051, 1057, 1059 } ;
     std::vector< uint64_t > L =  { 990, 396, 996, 60, 464, 468, 1008, 1012, 1018, 1020, 1030, 1032, 240, 1038,         346, 444, 1048, 1050, 150, 352  };
     std::vector< uint64_t > AB = {  126135, 125628, 125376, 124875, 124626, 124131, 123885, 123395, 122669,         122428, 121241, 121006, 120540, 120307, 120076, 119846, 119161, 118934, 118259,         118035   } ;
     std::vector< uint64_t > ub = { 10030183, 10016724, 10010022, 9996670, 9990021, 9976776, 9970180, 9957040, 9937459, 9930966, 9898754, 9892361, 9879625, 9873282, 9866955, 9860644, 9841808, 9835562, 9816916, 9810732 };
     
     std::vector< uint64_t > P =  {  1061, 1063, 1067, 1069, 1073, 1077, 1079, 1087, 1091, 1093, 1097, 1099, 1103,         1105, 1109, 1111, 1115, 1117, 1121, 1123 } ;
     std::vector< uint64_t > L =  { 1060, 1062, 480, 1068, 252, 358, 492, 1086, 1090, 1092, 1096, 156, 1102, 48,         1108, 100, 444, 1116, 522, 1122  };
     std::vector< uint64_t > AB = {  117813, 117591, 117150, 116931, 116495, 116063, 115848, 114995, 114573,         114364, 113947, 113739, 113327, 113122, 112714, 112511, 112107, 111906, 111507,         111308  } ;
     std::vector< uint64_t > ub = { 9804564, 9798411, 9786151, 9780045, 9767876, 9755769, 9749737, 9725760, 9713859, 9707931, 9696117, 9690232, 9678504, 9672661, 9661018, 9655217, 9643658, 9637898, 9626421, 9620703 };
         
     */
    
    for( int i = 0; i < 20; i++ )
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

old bigdawg 1 < P < 60      started 3-25:
new bigdawg 60 < P < 114    started 3-26: less than 48 hours
xeon phi    P = 115         started 3-27: 20 hours
xeon phi    P = 119         started 3-28: 22 hours
new bigdawg 120 < P < 178   started 8:21am 3-29: approx 24 hours
xeon phi    178 < P < 186   started 8:24am 3-29: approx 61 hours
nbd node-4  P = 187         started 2:27pm 3-30: approx 8 hours
xeon phi    P = 191         started 8:08am 4-1: 21 hours
nbd         192 < P < 240   started 7:40am 3-31: approx 25 hours
nbd         240 < P < 304   started 2:00pm 4-1: approx 27 hours
xeon phi    P = 307         started 9:25am 4-2: 17 hours
xeon phi    310 < P < 320   started 7:45am 4-3:
nbd         320 < P < 368   started 7:45am 4-3: about 16 hours
nbd node-0  370 < P < 427   started 6:00am 4-4:

magma script generated the inputs:

P:=[];  L:=[];  AB:= [];  ub:= [];  count:= 0;  i:= 1187;
while count lt 20 do
    if GCD( i, EulerPhi( i ) ) eq 1 then
        i_f:= Factorization( i ) ;
        if i*i_f[ #i_f ][1]^3 lt 10^24 then
            count:= count + 1;
            P[count]:= i;
            L[count]:= CarmichaelLambda(i);
            AB[count]:= Floor(125000000/i);
            ub[count]:= Ceiling( (10^24/i)^(1/3) + 1);
        end if;
    end if;
    i:= i + 2;
end while;
P;  L;  AB;  ub;
 
std::vector< uint64_t > P =   {        3,        5,        7,       11,      13,        15,       17,       19,       23,       29,       31,       33,       35,       37,       41,       43,       47,       51,       53,       59 };
std::vector< uint64_t > L =   {        2,        4,        6,       10,      12,         4,       16,       18,       22,       28,       30,       10,       12,       36,       40,       42,       46,       16,       52,       58 };
std::vector< uint64_t > AB =  { 41666666, 25000000, 17857142, 11363636,  9615384,  8333333,  7352941,  6578947,  5434782,  4310344,  4032258,  3787878,  3571428,  3378378,  3048780,  2906976,  2659574,  2450980,  2358490,  2118644 };
std::vector< uint64_t > ub =  { 69336128, 58480355, 52275796, 44964432, 42529038, 40548014, 38891112, 37475618, 35163390, 32548728, 31833138, 31176596, 30571072, 30010008, 29000494, 28543718, 27709842, 26965592, 26622044, 25687155 };
 

 std::vector< uint64_t > P =  {      371,      373,      377,      379,      383,      389,      391,      393,      395,      397,      401,      403,      407,      409,      411,      413,      415,      419,      421,      427 };
 std::vector< uint64_t > L =  {      156,      372,       84,      378,      382,      388,      176,      130,      156,      396,      400,       60,      180,      408,      136,      174,      164,      418,      420,       60 };
 std::vector< uint64_t > AB = {   336927,   335120,   331564,   329815,   326370,   321336,   319693,   318066,   316455,   314861,   311720,   310173,   307125,   305623,   304136,   302663,   301204,   298329,   296912,   292740 };
 std::vector< uint64_t > ub = { 13916886, 13891968, 13842661, 13818269, 13769995, 13698831, 13675434, 13652197, 13629116, 13606190, 13560798, 13538328, 13493830, 13471800, 13449912, 13428166, 13406560, 13363761, 13342566, 13279776 };
 
 std::vector< uint64_t > P =  {      431,      433,      435,      437,      439,      443,      445,      447,      449,      451,      455,      457,      461,      463,      467,      469,      473,      479,      481,      485 };
 std::vector< uint64_t > L =  {      430,      432,       28,      198,      438,      442,       88,      148,      448,       40,       12,      456,      460,      462,      466,       66,      210,      478,       36,       96 };
 std::vector< uint64_t > AB = {   290023,   288683,   287356,   286041,   284738,   282167,   280898,   279642,   278396,   277161,   274725,   273522,   271149,   269978,   267665,   266524,   264270,   260960,   259875,   257731 };
 std::vector< uint64_t > ub = { 13238566, 13218152, 13197863, 13177699, 13157656, 13117935, 13098253, 13078689, 13059241, 13039908, 13001584, 12982589, 12944931, 12926265, 12889253, 12870906, 12834521, 12780707, 12762968, 12727784 };

 std::vector< uint64_t > P =  {      487,      491,      493,      499,      501,      503,      509,      511,      515,      517,      519,      521,      523,      527,      533,     535,       537,      541,      545,      547 };
 std::vector< uint64_t > L =  {      486,      490,      112,      498,      166,      502,      508,       72,      204,      230,      172,      520,      522,      240,      120,     212,       178,      540,      108,      546 };
 std::vector< uint64_t > AB = {   256673,   254582,   253549,   250501,   249500,   248508,   245579,   244618,   242718,   241779,   240847,   239923,   239005,   237191,   234521,   233644,   232774,   231053,   229357,   228519 };
 std::vector< uint64_t > ub = { 12710337, 12675727, 12658563, 12607623, 12590824, 12574114, 12524511, 12508150, 12475682, 12459574, 12443549, 12427606, 12411744, 12380262, 12333632, 12318244, 12302932, 12272536, 12242437, 12227498 };
 
 std::vector< uint64_t > P =  {      551,      553,      557,      559,      561,      563,      565,      569,      571,      573,      577,      581,      583,      587,      589,      591,      593,      595,      599,      601 };
 std::vector< uint64_t > L =  {      252,       78,      556,       84,       80,      562,      112,      568,      570,      190,      576,      246,      260,      586,       90,      196,      592,       48,      598,      600 };
 std::vector< uint64_t > AB = {   226860,   226039,   224416,   223613,   222816,   222024,   221238,   219683,   218914,   218150,   216637,   215146,   214408,   212947,   212224,   211505,   210792,   210084,   208681,   207986 };
 std::vector< uint64_t > ub = { 12197838, 12183115, 12153881, 12139369, 12124926, 12110552, 12096245, 12067833, 12053727, 12039687, 12011801, 11984171, 11970452, 11943200, 11929666, 11916194, 11902782, 11889431, 11862907, 11849733 };
 
 */
