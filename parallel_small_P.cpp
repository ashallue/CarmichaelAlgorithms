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

    std::string file_name = "nbd_node1_1013_small_" + std::to_string(my_rank) + ".txt";

    uint16_t work_count = 0;

    std::vector< uint64_t > P =  {  991, 995, 997, 1001, 1003, 1007, 1009, 1013, 1019, 1021, 1031, 1033, 1037,        1039, 1041, 1043, 1049, 1051, 1057, 1059  };
    std::vector< uint64_t > L =  {   990, 396, 996, 60, 464, 468, 1008, 1012, 1018, 1020, 1030, 1032, 240, 1038,
        346, 444, 1048, 1050, 150, 352  };
    std::vector< uint64_t > AB = {  126135, 125628, 125376, 124875, 124626, 124131, 123885, 123395, 122669,
        122428, 121241, 121006, 120540, 120307, 120076, 119846, 119161, 118934, 118259,
        118035 };
    std::vector< uint64_t > ub = { 10030183, 10016724, 10010022, 9996670, 9990021, 9976776, 9970180, 9957040,
        9937459, 9930966, 9898754, 9892361, 9879625, 9873282, 9866955, 9860644, 9841808,
        9835562, 9816916, 9810732 };

    
    /*
    std::vector< uint64_t > P =  {    1495,    1499,    1501,    1507,    1509,    1511,    1513,    1517,    1523,    1527,    1529,    1531,    1535,    1537,    1541,    1543,    1547,    1549,    1551,    1553 };
    std::vector< uint64_t > L =  {     132,    1498,     234,     680,     502,    1510,     176,     360,    1522,     508,     690,    1530,     612,     364,      66,    1542,      48,    1548,     230,    1552 };
    std::vector< uint64_t > AB = {   83612,   83388,   83277,   82946,   82836,   82726,   82617,   82399,   82074,   81859,   81752,   81645,   81433,   81327,   81116,   81011,   80801,   80697,   80593,   80489 };
    std::vector< uint64_t > ub = { 8745534, 8737748, 8733866, 8722259, 8718404, 8714556, 8710714, 8703051, 8691608, 8684012, 8680224, 8676442, 8668899, 8665137, 8657633, 8653891, 8646426, 8642703, 8638987, 8635277 };
 
    std::vector< uint64_t > P =  {    1559,    1561,    1563,    1565,    1567,    1571,    1577,    1579,    1583,    1585,    1589,    1591,    1597,    1601,    1603,    1605,    1607,    1609,    1613,    1615 };
    std::vector< uint64_t > L =  {    1558,     222,     520,     312,    1566,    1570,     738,    1578,    1582,     316,     678,     252,    1596,    1600,     228,     212,    1606,    1608,    1612,     144 };
    std::vector< uint64_t > AB = {   80179,   80076,   79974,   79872,   79770,   79567,   79264,   79164,   78963,   78864,   78665,   78566,   78271,   78076,   77978,   77881,   77784,   77688,   77495,   77399 };
    std::vector< uint64_t > ub = { 8624184, 8620500, 8616821, 8613149, 8609483, 8602170, 8591246, 8587618, 8580378, 8576768, 8569565, 8565973, 8555232, 8548101, 8544544, 8540994, 8537449, 8533910, 8526850, 8523329 };
     
    std::vector< uint64_t > P =  {    1619,    1621,    1627,    1631,    1633,    1637,    1639,    1643,    1645,    1649,    1651,    1657,    1661,    1663,    1667,    1669,    1671,    1679,   1685,     1687 };
    std::vector< uint64_t > L =  {    1618,    1620,    1626,     696,     770,    1636,     740,     780,     276,      96,     252,    1656,     150,    1662,    1666,    1668,     556,     792,     336,     240 };
    std::vector< uint64_t > AB = {   77208,   77112,   76828,   76640,   76546,   76359,   76266,   76080,   75987,   75803,   75711,   75437,   75255,   75165,   74985,   74895,   74805,   74449,   74183,   74096 };
    std::vector< uint64_t > ub = { 8516303, 8512800, 8502322, 8495366, 8491896, 8484974, 8481521, 8474633, 8471197, 8464342, 8460923, 8450698, 8443909, 8440522, 8433766, 8430396, 8427031, 8413626, 8403627, 8400305 };

    std::vector< uint64_t > P =  {    1689,    1691,    1693,    1695,    1697,    1699,    1707,    1709,    1717,    1721,    1723,    1727,    1729,    1733,    1735,    1739,    1741,    1745,    1747,    1749 };
    std::vector< uint64_t > L =  {     562,     792,    1692,     112,    1696,    1698,     568,    1708,     400,    1720,    1722,     780,      36,    1732,     692,     828,    1740,     348,    1746,     260 };
    std::vector< uint64_t > AB = {   74008,   73920,   73833,   73746,   73659,   73572,   73227,   73142,   72801,   72632,   72547,   72379,   72296,   72129,   72046,   71880,   71797,   71633,   71551,   71469 };
    std::vector< uint64_t > ub = { 8396988, 8393676, 8390370, 8387068, 8383772, 8380481, 8367369, 8364103, 8351093, 8344618, 8341388, 8334943, 8331728, 8325313, 8322113, 8315727, 8312541, 8306185, 8303014, 8299848 };
     
     */
    
    for( int i = 7; i < 8; i++ )
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

old bigdawg     1 < P < 60          started 3-25: ended April 7th
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
obd             858 < P < 920       started 7:00am 4-8: about 2.5 days
nbd node-0      920 < P < 990       started 7:00am 4-8: a little over 2 days
nbd node-1      990 < P < 1060      started 12:10pm 4-9: about 42 hours - except 1 process (run again, see below)
nbd node-0      1060 < P < 1124     started 7:00am 4-10: about 41 hours
nbd node-2      1124 < P < 1184     started 7:00am 4-10: about 39 hours - except 1 process (!!! - run again)
nbd node-3      1184 < P < 1246     started 7:00am 4-10: about 38 hours
obd             1246 < P < 1362     started 8:00am 4-11:
nbd node-0      1362 < P < 1430     started 8:00am 4-12: less than 36 hours
nbd node-3      1430 < P < 1494     started 8:00am 4-12: less than 36 hours
xeon phi        P = 991             started 9:25am 4-14:
nbd node-0      P = 995             started 10:35am 4-14: done within 3 hours - error not here, all threads finished
nbd node-1      P = 997             started 11:10am 4-14: done within 3 hours - error not here, all threads finished
nbd node-2      P = 1001            started 12:15am 4-14:
nbd node-3      P = 1003            started 12:15am 4-14:
nbd node-4      P = 1007            started 12:15am 4-14:
nbd node-0      P = 1009            started 1:00pm, 4-14:
nbd node-1      P = 1013            started 1:10pm, 4-14:
 
 
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
 
std::vector< uint64_t > P =  {      991,      995,      997,    1001,    1003,    1007,    1009,    1013,    1019,    1021, 1031, 1033, 1037,        1039, 1041, 1043, 1049, 1051, 1057, 1059  };
std::vector< uint64_t > L =  {      990,      396,      996,      60,     464,     468,    1008,    1012,    1018,    1020, 1030, 1032, 240, 1038,     346, 444, 1048, 1050, 150, 352  };
std::vector< uint64_t > AB = {   126135,   125628,   125376,  124875,  124626,  124131,  123885,  123395,  122669,  122428, 121241, 121006, 120540, 120307, 120076, 119846, 119161, 118934, 118259,     118035 };
std::vector< uint64_t > ub = { 10030183, 10016724, 10010022, 9996670, 9990021, 9976776, 9970180, 9957040, 9937459, 9930966, 9898754, 9892361, 9879625, 9873282, 9866955, 9860644, 9841808, 9835562, 9816916, 9810732 };
  
 std::vector< uint64_t > P =  {    1247,    1249,    1253,    1257,    1259,    1261,    1267,    1271,    1273,    1277,    1279,    1283,    1285,    1289,    1291,    1293,    1295,    1297,    1301,    1303,
                                   1307,    1309,    1313,    1315,    1319,    1321,    *1327,    1329,    1333,    1335,    1337,    1339,    1343,    1345,    1347,    1349,    1351,    1353,    1357,    1361 };
 std::vector< uint64_t > L =  {      84,    1248,     534,     418,    1258,      96,     180,     120,     198,    1276,    1278,    1282,     256,    1288,    1290,     430,      36,    1296,    1300,    1302,
                                   1306,     240,     300,     524,    1318,    1320,    1326,     442,     210,      88,     570,     204,     624,     268,     448,     630,     192,      40,     638,    1360 };
 std::vector< uint64_t > AB = {  100240,  100080,   99760,   99443,   99285,   99127,   98658,   98347,   98193,   97885,   97732,   97427,   97276,   96974,   96824,   96674,   96525,   96376,   96079,   95932,
                                  95638,   95492,   95201,   95057,   94768,   94625,   94197,   94055,   93773,   93632,   93492,   93353,   93075,   92936,   92798,   92661,   92524,   92387,   92114,   91844 };
 std::vector< uint64_t > ub = { 9290618, 9285656, 9275764, 9265915, 9261006, 9256107, 9241473, 9231768, 9226931, 9217287, 9212480, 9202896, 9198119, 9188595, 9183847, 9179110, 9174382, 9169664, 9160257, 9155567,
                                9146218, 9141557, 9132265, 9127633, 9118396, 9113792, 9100036, 9095469, 9086362, 9081822, 9077291, 9072770, 9063753, 9059258, 9054772, 9050295, 9045827, 9041368, 9032475, 9023618 };
 
 std::vector< uint64_t > P =  {    1363,    1367,    1373,    1381,    1383,    1385,    1387,    1391,    1393,    1397,    1399,    1401,    1403,    1409,    1411,    1415,    1417,    1423,    1427,    1429 };
 std::vector< uint64_t > L =  {     644,    1366,    1372,    1380,     460,     276,      72,     636,     198,     630,    1398,     466,     660,    1408,     656,     564,     108,    1422,    1426,    1428 };
 std::vector< uint64_t > AB = {   91709,   91441,   91041,   90514,   90383,   90252,   90122,   89863,   89734,   89477,   89349,   89221,   89094,   88715,   88589,   88339,   88214,   87842,   87596,   87473 };
 std::vector< uint64_t > ub = { 9019202, 9010396, 8997252, 8979845, 8975514, 8971192, 8966878, 8958274, 8953985, 8945431, 8941166, 8936910, 8932661, 8919963, 8915747, 8907338, 8903145, 8890614, 8882300, 8878154 };

 
 */
