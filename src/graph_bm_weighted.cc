/*
 * ============================================================================
 *
 *       Filename:  graph_bm.cc
 *
 *         Author:  Prashant Pandey (), ppandey@lbl.gov
 *   Organization:  Berkeley Lab
 *
 * ============================================================================
 */

#define ENABLE_LOCK 1
#define WEIGHTED 1
#define VERIFY 0

#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <queue>
#include "sys/time.h"

#include "graph.h"
#include "parallel.h"
#include "util.h"
#include "parse_command_line.h"
#include "io_util.h"

#include "gabps/bitmap.h"
#include "gabps/platform_atomics.h"
#include "gabps/pvector.h"
// #include "gabps/sliding_queue.h"

#include "PMA.hpp"

#include "BellmanFord.h"
// #include "../include/DeltaStepping.h"

using namespace graphstore;

#define LOGGING_TICK (1ULL << 24)

#include	<stdlib.h>


std::string test_name[] = {
    "SSSP_BF",
};

/*
// from gap: random vals uniform in [1, 255]
#define MAXVAL 254
uint32_t rand_in_range(uint32_t max) { return rand() % max; }
*/
template <class Graph>
double execute(Graph& G, commandLine& P, std::string testname, int i) {
  if (testname == "SSSP_BF") {
    return test_sssp_bf(G, P);
  } else if (testname == "SSSP_DS") {
    return test_sssp_ds(G, P);
  } else {
    std::cout << "Unknown test: " << testname << ". Quitting." << std::endl;
    exit(0);
  }
}

template <class G>
double test_sssp_bf(G& GA, commandLine& P) {
	struct timeval start, end;
	struct timezone tzp;

  long src = P.getOptionLongValue("-src",-1);
  if (src == -1) {
    std::cout << "Please specify a source vertex to run the SSSP from" << std::endl;
    exit(0);
  }

  // std::cout << "Running SSSP BF from source = " << src << std::endl;
  // with edge map
  gettimeofday(&start, &tzp);
  auto result = SSSP_BF(GA, src);
  // SSSP_BF(GA, src);
  gettimeofday(&end, &tzp);
  // write out result to file for verification
  /*
  FILE *fp;
  fp = fopen("sssp.out", "w+");
  int32_t min_len = INT_MAX;
  for (int i = 0; i < GA.get_num_vertices(); i++) {
    fprintf(fp, "%d\n", result[i]);
  }
  fclose(fp);
  */

  free(result);

  return cal_time_elapsed(&start, &end);
}

template <class G>
double test_sssp_ds(G& GA, commandLine& P) {
	struct timeval start, end;
	struct timezone tzp;

  long src = P.getOptionLongValue("-src",-1);
  if (src == -1) {
    std::cout << "Please specify a source vertex to run the SSSP from" << std::endl;
    exit(0);
  }

  std::cout << "Running SSSP DS from source = " << src << std::endl;
  // with edge map
  gettimeofday(&start, &tzp);
  // auto result = SSSP_DS(GA, P); //src);
  gettimeofday(&end, &tzp);
  printf("SSSP DS finished\n");
  return cal_time_elapsed(&start, &end);
}

void run_algorithm(commandLine& P) {
  auto test_id = P.getOptionValue("-t", "SSSP_BF");
  size_t rounds = P.getOptionLongValue("-rounds", 4);
  double total_time = 0.0;
  
	std::string src, dest;
	// read edges as source and destination
	//int cnt = 0;
	struct timeval start, end;
	struct timezone tzp;

  // initialize graph
	uint32_t num_nodes;
  uint64_t num_edges;
  auto filename = P.getOptionValue("-f", "none");
	// trip_uint *edges = get_edges_from_file(filename.c_str(), 1,true, &num_edges, &num_nodes);
  trip_uint *edges = get_wgh_edges_from_file_adj_sym(filename.c_str(),&num_edges, &num_nodes);
	Graph graph(num_nodes);

//	std::random_device rd;
//	std::mt19937 g(rd());
//	std::shuffle(edges, edges+num_edges, g);

	std::vector<uint32_t> new_srcs(num_edges);
      std::vector<uint32_t> new_dests(num_edges);
      std::vector<uint32_t> new_wghts(num_edges);
             for (uint32_t i = 0; i < num_edges; i++) {
               new_srcs[i] = edges[i].x;
               new_dests[i] = edges[i].y;
							 new_wghts[i] = edges[i].z;
       }
       auto perm = get_random_permutation(num_edges);

	PRINT("Inserting edges");	

	gettimeofday(&start, &tzp);
  // randomly generate weights
	graph.add_edge_batch(new_srcs.data(), new_dests.data(), new_wghts.data(), num_edges, perm);
	//parallel_for_256 (uint64_t i = 0; i < num_edges; i++) {
	//for (uint64_t i = 0; i < num_edges; i++) {
    // uint32_t val = rand_in_range(MAXVAL) + 1;
		//graph.add_edge(edges[i].x, edges[i].y, edges[i].z);
    // fprintf(file, "%u\t%u\t%u\n", edges[i].x, edges[i].y, val);
  //}
	gettimeofday(&end, &tzp);
	free(edges);
	new_srcs.clear();
  new_dests.clear();
  new_wghts.clear();

	float size_gb = graph.get_size()/(float)1073741824;
	PRINT("Construction finished. Nodes: " << graph.get_num_vertices() <<
				" Edges: " << graph.get_num_edges() << " Size: " << size_gb << " GB");
	print_time_elapsed("Inserting edges: ", &start, &end);
	PRINT("Throughput: " <<
				graph.get_num_edges()/(float)cal_time_elapsed(&start, &end));

  for (size_t i=0; i<rounds; i++) {
    double tm = execute(graph, P, test_id, i);

    // std::cout << "RESULT"  << fixed << setprecision(6)
    std::cout << "\ttest=" << test_id
     << "\ttime=" << tm
     << "\titeration=" << i << std::endl;
    total_time += tm;
  }
  // std::cout << "RESULT (AVG)" << fixed << setprecision(6)
  std::cout << "AVG"
    << "\ttest=" << test_id
    << "\ttime=" << (total_time / rounds)
    << "\tgraph=" << filename << std::endl;
}

// take in unweighted graph, write out (symmetrized) weighted one
void output_weighted(commandLine& P) {
  auto filename = P.getOptionValue("-f", "none");

  // split
  std::vector<std::string> parts;
  std::stringstream ss(filename);
  std::string token;
  char delim = '/';
  while (std::getline(ss, token, delim)) {
      parts.push_back(token);
  }

  std::string outfile = parts[parts.size() - 1] + ".weighted";
	std::cout << "Writing weighted graph to " << outfile << std::endl;	

  auto edges = get_unique_edges_from_file(filename.c_str());

  FILE* file;
  file = fopen(outfile.c_str(), "w");
  uint64_t num_written = 0;
  // randomly generate weights
	for (auto itr = edges.begin(); itr != edges.end(); ++itr) {
    num_written++;
    fprintf(file, "%u\t%u\t%u\n",itr->first.first, itr->first.second, itr->second);
    fprintf(file, "%u\t%u\t%u\n",itr->first.second, itr->first.first, itr->second);
  }
  printf("num written %lu\n", num_written);
  fclose(file);
}

/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
int main(int argc, char** argv) {
  srand(time(NULL));
  printf("Num workers: %ld\n", getWorkers());
  printf("Num workers: %ld\n", getWorkers());
  commandLine P(argc, argv, "./graph_bm [-t testname -r rounds -f file");
  auto todo = P.getOptionValue("-t", "none");

  // take in an unweighted graph and weight it, write it to ../inputs
  if (todo == "WEIGHT") {
    output_weighted(P);
  } else {
    std::cout <<  todo << std::endl;
    // otherwise, run algorithm on weighted graph
    run_algorithm(P);
  }
}
