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
#define WEIGHTED 0
#define VERIFY 0

#include <stdlib.h>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <queue>
#include "sys/time.h"

#include "graph.h"
#include "parallel.h"
#include "util.h"
#include "parse_command_line.h"
#include "rmat_util.h"

#include "gabps/bitmap.h"
#include "gabps/platform_atomics.h"
#include "gabps/pvector.h"
// #include "gabps/sliding_queue.h"

#include "integerSort/blockRadixSort/blockRadixSort.h"
#include "PMA.hpp"

#include "BFS.h"
#include "Pagerank.h"
#include "Components.h"
#include "BC.h"
#include "TC.h"
#include "BellmanFordUnweighted.h"
#include "io_util.h"

using namespace graphstore;

#define LOGGING_TICK (1ULL << 24)
#define BATCH_SIZE (1ULL << 9)

#include	<stdlib.h>

std::string test_name[] = {
    "BFS",
    "PR",
    "CC",
    "BC",
    "SSSP_BF",
};

template <class Graph>
double execute(Graph& G, commandLine& P, std::string testname, int i) {
  if (testname == "BFS") {
    return test_bfs(G, P,i );
  } else if (testname == "PR") {
    return test_pr(G, P);
  } else if (testname == "CC") {
    return test_cc(G, P);
  } else if (testname == "BC") {
    return test_bc(G, P);
  } else if (testname == "TC") {
    return test_tc(G, P);
  } else if (testname == "SSSP_BF") {
    return test_sssp_bf(G, P);
  } else {
    std::cout << "Unknown test: " << testname << ". Quitting." << std::endl;
    exit(0);
  }
}

template <class G>
double test_pr(G& GA, commandLine& P) {
	struct timeval start, end;
	struct timezone tzp;

  long maxiters = P.getOptionLongValue("-maxiters",10);

  std::cout << "Running PR" << std::endl;

  // with edge map
  gettimeofday(&start, &tzp);
  auto pr_edge_map = PR_S<double>(GA, maxiters); 
  gettimeofday(&end, &tzp);
  free(pr_edge_map);
  PRINT("PR finished");
  return cal_time_elapsed(&start, &end);
}

template <class G>
double test_cc(G& GA, commandLine& P) {
	struct timeval start, end;
	struct timezone tzp;

  std::cout << "Running CC" << std::endl;

  // with edge map
  gettimeofday(&start, &tzp);
  auto cc_result = CC(GA);

  /*
  std::ofstream myfile;
  myfile.open ("cc.out");
  for (int i = 0; i < GA.get_num_vertices(); i++) {
    myfile << cc_result[i] << "\n";
  }
  myfile.close();
  */
  gettimeofday(&end, &tzp);
  free(cc_result);
  PRINT("CC finished");
  return cal_time_elapsed(&start, &end);
}

template <class G>
double test_bc(G& GA, commandLine& P) {
	struct timeval start, end;
	struct timezone tzp;

  long src = P.getOptionLongValue("-src",-1);
  if (src == -1) {
    std::cout << "Please specify a source vertex to run the BC from" << std::endl;
    exit(0);
  }
  std::cout << "Running BC from source = " << src << std::endl;

  // with edge map
  gettimeofday(&start, &tzp);
  auto bc_result = BC(GA, src);
  gettimeofday(&end, &tzp);
  free(bc_result);
  PRINT("BC finished");
  return cal_time_elapsed(&start, &end);
}

template <class G>
double test_tc(G& GA, commandLine& P) {
	struct timeval start, end;
	struct timezone tzp;

  std::cout << "Running TC" << std::endl;

  // with edge map
  gettimeofday(&start, &tzp);
  auto count = TC(GA);
  gettimeofday(&end, &tzp);
  printf("TC finished, counted %ld\n", count);
  return cal_time_elapsed(&start, &end);
}


// return time elapsed
template <class G>
double test_bfs(G& GA, commandLine& P, int trial) {
	struct timeval start, end;
	struct timezone tzp;

  long src = P.getOptionLongValue("-src",-1);
  if (src == -1) {
    std::cout << "Please specify a source vertex to run the BFS from" << std::endl;
    exit(0);
  }
  std::cout << "Running BFS from source = " << src << std::endl;

  // without edge map
  /*
  auto num_edges = GA.get_num_edges();
  gettimeofday(&start, &tzp);
  pvector<int32_t> bfs_result = parallel_bfs(GA, src, num_edges);
  PRINT("BFS finished. parent of " << src << ": " << bfs_result[src]);
  gettimeofday(&end, &tzp);
  print_time_elapsed("parallel BFS: ", &start, &end);
  */

  // with edge map
  gettimeofday(&start, &tzp);
  auto bfs_edge_map = BFS_with_edge_map(GA, src);
  gettimeofday(&end, &tzp);
  PRINT("BFS finished. parent of " << src << ": " << bfs_edge_map[src]);
  free(bfs_edge_map);

  return cal_time_elapsed(&start, &end);
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
  gettimeofday(&end, &tzp);
  /*
  FILE *fp;
  fp = fopen("test.out", "w+");
  int32_t min_len = INT_MAX;
  for (int i = 0; i < GA.get_num_vertices(); i++) {
    min_len = std::min(min_len, result[i]);
    fprintf(fp, "%d\n", result[i]);
  }
  fclose(fp);
  */
  free(result);

  return cal_time_elapsed(&start, &end);
}


void run_algorithm(commandLine& P) {
  auto test_id = P.getOptionValue("-t", "BFS");
  size_t rounds = P.getOptionLongValue("-rounds", 4);
  
	std::string src, dest;
	// read edges as source and destination
	//int cnt = 0;
	struct timeval start, end;
	struct timezone tzp;

  // initialize graph
	uint32_t num_nodes;
  uint64_t num_edges;
  auto filename = P.getOptionValue("-f", "none");
	//pair_uint *edges = get_edges_from_file(filename.c_str(), 1,true, &num_edges, &num_nodes);
	pair_uint *edges = get_edges_from_file_adj_sym(filename.c_str(), &num_edges, &num_nodes);

	Graph graph(num_nodes);
	//std::random_shuffle(edgelist.begin(), edgelist.end());

	//std::random_device rd;
	//std::mt19937 g(rd());
	//std::shuffle(edges, edges+num_edges, g);
      std::vector<uint32_t> new_srcs(num_edges);
      std::vector<uint32_t> new_dests(num_edges);
             for (uint32_t i = 0; i < num_edges; i++) {
               new_srcs[i] = edges[i].x;
               new_dests[i] = edges[i].y;
       }
       auto perm = get_random_permutation(num_edges);

	PRINT("Inserting edges");	
	gettimeofday(&start, &tzp);
	graph.add_edge_batch(new_srcs.data(), new_dests.data(), num_edges, perm);
	//parallel_for (uint64_t i = 0; i < num_edges; i++) {
	//for (uint64_t i = 0; i < num_edges; i++) {
		//graph.add_edge(edges[i].x, edges[i].y);
	//}

	gettimeofday(&end, &tzp);
  for(uint32_t i = 0; i < num_edges; i++) {
    if (!graph.is_edge(edges[i].x, edges[i].y)) {
      printf("edge (%u, %u) not found, should be\n", edges[i].x, edges[i].y);
      printf("\tdegree %u = %u\n", edges[i].x, graph.degree(edges[i].x));
    }
  }

	free(edges);
	       new_srcs.clear();
       new_dests.clear();
	float size_gb = graph.get_size()/(float)1073741824;
	PRINT("Construction finished. Nodes: " << graph.get_num_vertices() <<
				" Edges: " << graph.get_num_edges() << " Size: " << size_gb << " GB");
	print_time_elapsed("Inserting edges: ", &start, &end);
	PRINT("Throughput: " <<
				graph.get_num_edges()/(float)cal_time_elapsed(&start, &end));
//#if 0
  std::vector<std::string> test_ids;
  // if testname is TC, include it, otherwise exclude it
  if (P.getOptionLongValue("-TC",0) != 0) {
    test_ids = {"BFS","PR","CC","BC", "TC"};
  } else {
    test_ids = {"BFS","PR","CC","BC"};
  }
  for (auto test_id : test_ids) {
    double total_time = 0.0;
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
  commandLine P(argc, argv, "./graph_bm [-t testname -r rounds -f file");
  run_algorithm(P);
}
/*
	int
main ( int argc, char *argv[] )
{
	srand(time(NULL));
#if 0
    uint32_t num_nodesPMA = 100;
	 PMA pma = PMA(num_nodesPMA);
	 for (int i = 0; i < 50; i++) {
		 for (int j = 0; j < 100; j++) {
        uint32_t y = rand()%num_nodesPMA;
        //printf("adding edge (%u, %u\n", i, y);
			 pma.add_edge(i,y,1);
       //pma.print_array();
       if (rand()%5 == 0) {
        //printf("removing edge (%u, %u\n", i, y);
        pma.remove_edge(i,y);
        //pma.print_array();
       }
		 }
	 }
   //return 0;


	 for (int i = 0; i < 10; i++) {
		 std::cout << "Neighbors for node: " << i << " ";
		 for (auto it = pma.begin(i); it != pma.end(i); ++it) {
			 std::cout << (*it).dest << " ";
		 }
		 std::cout << "\n";
	 }

	 //pma.print_array(0);

#else
	
	if (argc < 3) {
		fprintf(stderr, "Please specify the name of the input file and node id\n");
		exit(1);
	}

	std::string src, dest;
	// read edges as source and destination
	//int cnt = 0;
	struct timeval start, end;
	struct timezone tzp;

	uint32_t num_nodes;
  	uint64_t num_edges;
	pair_uint *edges = get_edges_from_file(argv[1], 1,true, &num_edges, &num_nodes);

	Graph graph(num_nodes);
	//std::random_shuffle(edgelist.begin(), edgelist.end());

	PRINT("Inserting edges");	
	gettimeofday(&start, &tzp);

	for (uint64_t i = 0; i < num_edges; i++) {
		graph.add_edge(edges[i].x, edges[i].y);
	}
	gettimeofday(&end, &tzp);
	free(edges);

	PRINT("Construction finished. Nodes: " << graph.get_num_vertices() <<
				" Edges: " << graph.get_num_edges());
	print_time_elapsed("Inserting edges: ", &start, &end);
	PRINT("Throughput: " <<
				graph.get_num_edges()/(float)cal_time_elapsed(&start, &end));

  PRINT("Starting BFS");
  gettimeofday(&start, &tzp);
  auto bfs_result = BFS_with_edge_map(graph, atoi(argv[2]));
  //pvector<int32_t> bfs_result = parallel_bfs(graph, atoi(argv[2]), num_edges);
  //std::vector<int32_t> bfs_result = bfs(graph, atoi(argv[2]));
  gettimeofday(&end, &tzp);
  PRINT("BFS finished. parent of 0: " << bfs_result[0]);
  print_time_elapsed("BFS traversal: ", &start, &end);
#if 0
	PRINT("Starting BFS");
	auto it = graph.find();
	uint64_t total = 0;
	gettimeofday(&start, &tzp);
	while (!it.done()) {
		total++;
		++it;
	}
	gettimeofday(&end, &tzp);
	PRINT("BFS finished. Nodes: " << total);
	print_time_elapsed("BFS traversal: ", &start, &end);
#endif

#endif
	return EXIT_SUCCESS;
}*/				/* ----------  end of function main  ---------- */
