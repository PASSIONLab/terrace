/*
 * ============================================================================
 *
 *       Filename:  test_graphcontainer.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
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
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <openssl/rand.h>
#include "sys/time.h"

#include "integerSort/blockRadixSort/blockRadixSort.h"
#include "util.h"
#include "io_util.h"
#include "spdlog/spdlog.h"
//#include "graph_qf.h"
#include "graph.h"

using namespace graphstore;
std::shared_ptr<spdlog::logger> console;

#if WEIGHTED
// three tuple uint
//typedef struct trip_uint {
  //uint32_t x;
  //uint32_t y;
  //uint32_t z;
//} trip_uint;

#if 0
trip_uint * get_edges_from_file(const char *filename, int zero_indexed, bool add_both_directions, uint64_t *edge_count, uint32_t *node_count) {
  printf("counting lines in file %s\n", filename);
  FILE *fp;
  fp = fopen(filename, "r");
  setvbuf ( fp , NULL , _IOFBF , 1<<20 );
  size_t buf_size = 64;
  char *line = (char *) malloc(buf_size);
  uint64_t line_count = 0;
  if (fp) {
    while (getline(&line, &buf_size, fp) != -1) {
      if (line[0] == '#') continue;
      line_count++;
    }
    // return 0;
  } else {
    printf("file was not opened\n");
    exit(EXIT_FAILURE);
  }
  if (add_both_directions) {
    line_count *= 2;
  }
  *edge_count = line_count;
  printf("getting edges from file %s, add both sides %d\n", filename, add_both_directions);
  rewind(fp);
  trip_uint *edges = (trip_uint *) malloc(line_count * sizeof(trip_uint));
  uint32_t num_nodes = 0;
  uint64_t index = 0;
  while (getline(&line, &buf_size, fp) != -1) {
    if (line[0] == '#') continue;
    uint32_t elem_1;
    uint32_t elem_2;
    uint32_t elem_3;
    sscanf(line, "%u   %u   %u", &elem_1, &elem_2, &elem_3);
    num_nodes = std::max(num_nodes, elem_1);
    num_nodes = std::max(num_nodes, elem_2);
    uint32_t src = elem_1 + zero_indexed;
    uint32_t dest = elem_2 + zero_indexed;
    num_nodes = std::max(num_nodes, src);
    num_nodes = std::max(num_nodes, dest);
    src--;
    dest--;
    edges[index++] = {src,dest, elem_3};
    if (add_both_directions) {
      edges[index++] = {dest, src, elem_3};
    }
  }
  printf("weighted num edges %lu\n", index);
  fclose(fp);
    // return 0;
  free(line);
  *node_count = num_nodes;
  return edges;
}
#endif
#endif

#if 0
pair_uint * get_edges_from_file(const char *filename, int zero_indexed, bool add_both_directions, uint64_t *edge_count, uint32_t *node_count) {
  printf("counting lines in file %s\n", filename);
  FILE *fp;
  fp = fopen(filename, "r");
  setvbuf ( fp , NULL , _IOFBF , 1<<20 );
  size_t buf_size = 64;
  char *line = (char *) malloc(buf_size);
  uint64_t line_count = 0;
  if (fp) {
    while (getline(&line, &buf_size, fp) != -1) {
      if (line[0] == '#') continue;
      line_count++;
    }
    // return 0;
  } else {
    printf("file was not opened\n");
    exit(EXIT_FAILURE);
  }
  if (add_both_directions) {
    line_count *= 2;
  }
  *edge_count = line_count;
  printf("getting edges from file %s\n", filename);
  rewind(fp);
  pair_uint *edges = (pair_uint *) malloc(line_count * sizeof(pair_uint));
  uint32_t num_nodes = 0;
  uint64_t index = 0;
  while (getline(&line, &buf_size, fp) != -1) {
    if (line[0] == '#') continue;
    uint32_t elem_1;
    uint32_t elem_2;
    sscanf(line, "%u   %u", &elem_1, &elem_2);
    num_nodes = std::max(num_nodes, elem_1);
    num_nodes = std::max(num_nodes, elem_2);
    uint32_t src = elem_1 + zero_indexed;
    uint32_t dest = elem_2 + zero_indexed;
    num_nodes = std::max(num_nodes, src);
    num_nodes = std::max(num_nodes, dest);
    src--;
    dest--;
    edges[index++] = {src,dest};
    if (add_both_directions) {
      edges[index++] = {dest, src};
    }
  }
  fclose(fp);
    // return 0;
  free(line);
  *node_count = num_nodes;
  return edges;
}
#endif

struct pair_hash
{
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2> &pair) const
	{
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

int test_correctness(Graph& graph, std::unordered_map<uint32_t,
                      Graph::vertex_set>& adj_list,
											std::unordered_map<std::pair<uint32_t, uint32_t>,
											uint32_t, pair_hash>& edge_map) {
  // check graph correctness by iterating over @adj_list
  for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
    //std::cout << "Degree: " << graph.degree(it->first) << "\n";
    if (graph.degree(it->first) != it->second.size()) {
      ERROR("Degree check failed for node: " << it->first << " Actual: " <<
            it->second.size() << " Reported: " << graph.degree(it->first));
      return EXIT_FAILURE;
    }
    auto neighbors = graph.neighbors(it->first);
    while (!neighbors.done()) {
      //std::cout << *neighbors << " ";
#if WEIGHTED
      auto ret = it->second.find((*neighbors).first);
#else
      auto ret = it->second.find(*neighbors);
#endif
      if (ret == it->second.end()) {
        ERROR("Neighbor check failed for node: " << it->first );
        return EXIT_FAILURE;
      }
      ++neighbors;
    }
    //PRINT("");
  }
  PRINT("Correctness check passed: by iterating over @adj_list!");

  // check graph correctness by iterating over @edge_map
  for (auto it = edge_map.begin(); it != edge_map.end(); ++it) {
#if WEIGHTED
    auto weight = graph.is_edge((*it).first.first, (*it).first.second);
    if (weight != (*it).second) {
      PRINT("Graph edge: " << (*it).first.first << " - " << (*it).first.second);
      PRINT("Weight reported: " << weight << " actual: " << (*it).second);
#else
    if (!graph.is_edge((*it).first.first, (*it).first.second)) {
      PRINT("Graph edge: " << (*it).first.first << " - " << (*it).first.second);
#endif
      ERROR("Edge test check failed!");
      return EXIT_FAILURE;
    }
  }
  PRINT("Correctness check passed: by iterating over @edge_list!");

  return EXIT_SUCCESS;
}

/* 
 * ===  FUNCTION  =============================================================
 *         Name:  main
 *  Description:  
 * ============================================================================
 */
  int
main ( int argc, char *argv[] )
{
  if (argc < 2) {
    fprintf(stderr, "Please specify the log of the number of slots in the CQF.\n");
    exit(1);
  }
#if WEIGHTED 
	PRINT("Testing weighted.");
#endif
  uint64_t qbits{0};
  uint64_t degree{0};
  uint64_t nslots{0};
  uint64_t nvals{0};
  uint32_t num_nodes{0};
  uint64_t num_edges{0};
  bool file = false;
#if WEIGHTED
    trip_uint *edges;
#else
  pair_uint *edges;
#endif
  if (argc == 2) {
    file = true;
#if WEIGHTED
    //edges = get_edges_from_file(argv[1], 1,true, &num_edges, &num_nodes);
  	edges = get_wgh_edges_from_file_adj_sym(argv[1],&num_edges, &num_nodes);
#else
    //edges = get_edges_from_file(argv[1], 1,true, &num_edges, &num_nodes);
		edges = get_edges_from_file_adj_sym(argv[1], &num_edges, &num_nodes);
#endif
  nslots = num_nodes;
  } else if (argc == 3) {
    qbits = atoi(argv[1]);
    degree = atoi(argv[2]);
    nslots = (1ULL << qbits);
		num_nodes = nslots;
    nvals = 100*nslots/100;
  }
	
	//std::random_device rd;
	//std::mt19937 g(rd());
	//std::shuffle(edges, edges+num_edges, g);

  //std::ofstream graph_file("data/random_graph.shuf", std::ofstream::out); 

  // create a typedef for the Graph type

  console = spdlog::default_logger();
  Graph graph(nslots);

	std::cout << "Inserting edges\n";
  std::unordered_map<Graph::vertex, Graph::vertex_set> adj_list;
	std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t, pair_hash> edge_map;
  //Graph::edge_set edge_list;
  if (!file) {
		/* Generate random values */
		uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
		//RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
		//for (uint32_t i = 0; i < nvals; i++) {
			//vals[i] = (1 * vals[i]) % nvals;
		//}

		srand(0);
		for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (rand() % nvals);
		}

		std::ofstream graph_fl("data/random_graph.weighted");
		std::ofstream graph_dot("data/random_graph.weighted.dot");
		graph_dot << "graph weighted {\n";
		for (uint32_t i = 0; i < nvals; i++) {
			uint32_t key = vals[i];
			// if not already seen
			if (adj_list.find(key) == adj_list.end()) {
				uint32_t nedges = rand() % degree + 1;

				Graph::vertex_set vec;
				for (uint32_t j = 0; j < nedges; j++) {
					uint32_t tonode = vals[rand() % nvals];
					if (adj_list[tonode].find(key) == adj_list[tonode].end() && key !=
							tonode) {
						if (vec.insert(tonode).second) {
#if WEIGHTED
							uint32_t w = vals[rand() % nvals] + 1;
							//graph.add_edge(key, tonode, w);
							edge_map[{key, tonode}] = w;
							graph_fl << key << "\t" << tonode << "\t" << w << '\n';
							graph_dot << key << " -- " << tonode << "[label=" << std::to_string(w) << "]\n";
#else
							//graph.add_edge(key, tonode);
							edge_map[{key, tonode}] = 1;
							graph_fl << key << "\t" << tonode << '\n';
							graph_dot << key << " -- " << tonode << "\n";
#endif
						}
					}
				}
				adj_list[key].merge(vec);
			}
		}
		graph_dot << '}';
		graph_fl.close();
		graph_dot.close();
  } else {
    for (uint32_t i = 0; i < num_edges; i++) {
			std::cout << "\r" << i << std::flush;
      auto src = edges[i].x;
      auto dest = edges[i].y;
#if WEIGHTED
      auto w = edges[i].z;
      //graph.add_edge(src, dest, w);
			edge_map[{src, dest}] = w;
#else
      //graph.add_edge(src, dest);
			edge_map[{src, dest}] = 1;
#endif
      adj_list[src].insert(dest);
    }
  }

	struct timeval start, end;
	struct timezone tzp;

	uint32_t *srcs = (uint32_t*)calloc(num_edges, sizeof(uint32_t));
	uint32_t *dests = (uint32_t*)calloc(num_edges, sizeof(uint32_t));
#if WEIGHTED
	uint32_t *wghts = (uint32_t*)calloc(num_edges, sizeof(uint32_t));
#endif
	if (!file) {
		uint64_t num_edges = edge_map.size();
		pair_uint *edges = (pair_uint*)calloc(num_edges, sizeof(pair_uint));
		uint32_t i = 0;
		for (auto edge : edge_map) { 
				edges[i].x = edge.first.first;
				edges[i].y = edge.first.second;
				i++;
		}
		integerSort_y((pair_els*)edges, num_edges, num_nodes);
		integerSort_x((pair_els*)edges, num_edges, num_nodes);
		for (uint32_t i = 0; i < num_edges; i++) {
			srcs[i] = edges[i].x;
			dests[i] = edges[i].y;
		}
		//for (auto edge : edge_map) 
		//edge_vector.push_back({{edge.first.first, edge.first.second}, edge.second});
		auto perm = get_random_permutation(num_edges);
		gettimeofday(&start, &tzp);
#if !WEIGHTED
		graph.add_edge_batch(srcs, dests, num_edges, perm);
#endif
		//for (uint32_t i = 0; i < edge_vector.size(); i++) {
			//auto edge = edge_vector[i];
//#if WEIGHTED
			//graph.add_edge(edge.first.first, edge.first.second, edge.second);
//#else
			//graph.add_edge(edge.first.first, edge.first.second);
//#endif
		//}
	} else {
		for (uint32_t i = 0; i < num_edges; i++) {
			srcs[i] = edges[i].x;
			dests[i] = edges[i].y;
#if WEIGHTED
			wghts[i] = edges[i].z;
#endif
		}
		auto perm = get_random_permutation(num_edges);
		gettimeofday(&start, &tzp);
#if WEIGHTED
		graph.add_edge_batch(srcs, dests, wghts, num_edges, perm);
#else
		graph.add_edge_batch(srcs, dests, num_edges, perm);
#endif
		//for (uint32_t i = 0; i < num_edges; i++) {
//#if WEIGHTED
			//graph.add_edge(edges[i].x, edges[i].y, edges[i].w);
//#else
			//graph.add_edge(edges[i].x, edges[i].y);
//#endif
		//}
	}
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);

  //graph_file.close();

  PRINT("Num vertices: " << graph.get_num_vertices());
  PRINT("Num edges: " << graph.get_num_edges());
  PRINT("Size in Bytes: " << graph.get_size());

	if (test_correctness(graph, adj_list, edge_map) == EXIT_FAILURE)
		return EXIT_FAILURE;

	std::cout << "Deleting edges\n";
	uint32_t num_prev{0};
	if (!file) {
		gettimeofday(&start, &tzp);
		for (uint32_t i = 0; i < num_edges; i++) {
			graph.remove_edge(srcs[i], dests[i]);
			if (num_prev == 0)
				num_prev = graph.get_num_edges();
			else if (num_prev - 1 != graph.get_num_edges()) {
				PRINT("Num edges: " << graph.get_num_edges() << " " << srcs[i] << " " << dests[i]);
				graph.print_vertex_block(srcs[i]);
				return EXIT_FAILURE;
			} else {
				num_prev = graph.get_num_edges();
			}
		}
	} else {
		gettimeofday(&start, &tzp);
		for (uint32_t i = 0; i < num_edges; i++) {
			graph.remove_edge(srcs[i], dests[i]);
		}
	}
	gettimeofday(&end, &tzp);
	print_time_elapsed("", &start, &end);

  PRINT("Num vertices: " << graph.get_num_vertices());
  PRINT("Num edges: " << graph.get_num_edges());
  PRINT("Size in Bytes: " << graph.get_size());

  return EXIT_SUCCESS;
}       /* ----------  end of function main  ---------- */
