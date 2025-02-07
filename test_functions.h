/*  Test functions for tabulating carmichaels project.
 *  Andrew Shallue and Jonathan Webster, January 2025
 */

#include "Preproduct.h"
#include <iostream>
#include <fstream>

// call CN_factorization on all squarefree multiples of 7 with 4 prime factors
bool test_factor();

// function that tabulates all Carmichaels up to a given bound.  Only large case, no small case
// need to read in jobs from output_jobs.txt, create a preproduct object for each one, call CN search, assemble the output
void tabulate_test(uint64_t bound, std::string jobs_file, std::string cars_file);
