/*  Test functions for tabulating carmichaels project.
 *  Andrew Shallue and Jonathan Webster, January 2025
 */

#include "Preproduct.h"
#include "IncrementalSieve/rollsieve.h"
#include <iostream>
#include <fstream>
#include <chrono>

// call CN_factorization on all squarefree multiples of 7 with 4 prime factors
bool test_factor();

// function that tabulates all Carmichaels up to a given bound.  Only large case, no small case
// need to read in jobs from output_jobs.txt, create a preproduct object for each one, call search, assemble the output
// Currently the method variable would be CN_search or CN_search_no_wheel
void tabulate_test(uint64_t bound, std::string jobs_file, std::string cars_file);

// function that prints wall time for one particular job
void job_timing(uint64_t P, uint64_t L, uint64_t prime_lower, std::string cars_file);

/*

*/
