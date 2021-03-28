/*
 * ============================================================================
 *
 *       Filename:  test_pma.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <random>
#include <algorithm>
#include <queue>
#include "sys/time.h"
#include <openssl/rand.h>

#include "util.h"

#include "PMA.hpp"

#include "util.h"
#include "io_util.h"

#include	<stdlib.h>

#define VERIFY 0

using namespace graphstore;

/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
	int
main ( int argc, char *argv[] )
{
	srand(time(NULL));
  uint32_t num_nodes{0};
  uint64_t num_edges{0};
  pair_uint *edges = get_edges_from_file_adj_sym(argv[1], &num_edges, &num_nodes);
	PMA pma = PMA(num_nodes);
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(edges, edges+num_edges, g);
	
	struct timeval start, end;
	struct timezone tzp;

	gettimeofday(&start, &tzp);
	for (uint32_t i=0; i<num_edges; i++) {
		pma.add_edge_update(edges[i].x, edges[i].y, 1);
	}
	gettimeofday(&end, &tzp);
	print_time_elapsed("Insertion time: ", &start, &end);

	gettimeofday(&start, &tzp);
	for (uint32_t i=0; i<num_edges; i++) {
		if (!pma.find_value(edges[i].x, edges[i].y))	{
			ERROR("Query failed!");
			return EXIT_FAILURE;
		}
	}
	gettimeofday(&end, &tzp);
	print_time_elapsed("Query time: ", &start, &end);
	
	gettimeofday(&start, &tzp);
	for (uint32_t i=0; i<num_edges; i++) {
		if (!pma.remove_edge(edges[i].x, edges[i].y))	{
			ERROR("Query failed!");
			return EXIT_FAILURE;
		}
	}
	gettimeofday(&end, &tzp);
	print_time_elapsed("Deletion time: ", &start, &end);

	PRINT("Test passed!!");
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
