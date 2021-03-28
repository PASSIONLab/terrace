/*
 * ============================================================================
 *
 *       Filename:  test_btree.cc
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#define WEIGHTED 1

#include <stdlib.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <openssl/rand.h>

#include "btree.h"
#include "util.h"

	struct Map {

#if WEIGHTED 
		std::set<std::pair<uint32_t, uint32_t>> set;
#else
		std::set<uint32_t> set;
#endif
		Map() {}

#if WEIGHTED
		void update(uint32_t v, uint32_t w) {
			set.insert({v,w});
		}
#else
		void update(uint32_t v) {
			set.insert(v);
		}
#endif
	};


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
		fprintf(stderr, "Please specify the log of the number of items.\n");
		exit(1);
	}
	uint64_t qbits = atoi(argv[1]);
	uint64_t nslots = (1ULL << qbits);
	uint64_t nvals = 100*nslots/100;

#if WEIGHTED 
	PRINT("Testing weighted.");
#endif

	srand(0);
	/* Generate random values */
	uint32_t *vals = (uint32_t*)malloc(nvals*sizeof(vals[0]));
	//RAND_bytes((unsigned char *)vals, sizeof(*vals) * nvals);
	//for (uint32_t i = 0; i < nvals; i++) {
		//vals[i] = (1 * vals[i]);
	//}
	for (uint32_t i = 0; i < nvals; i++) {
		vals[i] = (rand() % nvals);
		//vals[i] = i;
	}

	BTree<uint32_t, uint32_t> *tree = new BTree<uint32_t, uint32_t>();
#if WEIGHTED
	std::map<uint32_t, uint32_t> set;
#else
	std::set<uint32_t> set;
#endif

	for (uint32_t i = 0; i < nvals; i++) {
#if WEIGHTED
		auto it = set.insert(std::make_pair(vals[i], (uint32_t)(vals[i+1]%nvals)));
#else
		auto it = set.insert(vals[i]);
#endif

		if (it.second) {
#if WEIGHTED
			bool ret = tree->insert(vals[i], vals[i+1]%nvals);
#else
			bool ret = tree->insert(vals[i]);
#endif
			if (!ret) {
				ERROR("Failed return value for " << vals[i]);
#if WEIGHTED
				ERROR("Weight in Btree: " << tree->get_val(vals[i]));
				ERROR("Weight in set: " << (uint32_t)(vals[i+1]%nvals));
#endif
				ERROR("Btree: " << ret << " Set: " << it.second);
				//return EXIT_FAILURE;
			}
		}
		//std::cout << vals[i] << '\n';
		//for (auto it = tree->begin(); !it.done(); ++it) {
		//std::cout << *it << " ";
		//}
		//PRINT("");
	}

	//std::cout << "\nQuery:\n";
	for (uint32_t i = 0; i < nvals; i++) {
		//std::cout << vals[i] << ' ';
		if (tree->find(vals[i]) == nullptr) {
			ERROR("Query failed!");
			return EXIT_FAILURE;
		}
	}

	for (auto k : set) {
		//std::cout << k << " ";
#if WEIGHTED
		if (tree->find(k.first) == nullptr) {
#else
		if (tree->find(k) == nullptr) {
#endif
			ERROR("Query failed!");
			return EXIT_FAILURE;
		}
	}
	//PRINT("");

	//tree->traverse();
	//PRINT("");

	for (auto it = tree->begin(); !it.done(); ++it) {
#if WEIGHTED
		auto ret = set.find((*it).first);
#else
		auto ret = set.find(*it);
#endif
		//std::cout << *it << " ";
		if (ret == set.end()) {
			ERROR("Query failed!");
			return EXIT_FAILURE;
		}
	}
	//PRINT("");
	
#if !WEIGHTED
	Map map;
	tree->map(map);

	if (map.set != set) {
		ERROR("Query failed!");
		return EXIT_FAILURE;
	}
#endif

	//std::cout << "\nRemove\n";
	for (uint32_t i = 0; i < nvals; i++) {
		//std::cout << vals[i] << ' ';
		if (tree->find(vals[i]) != nullptr) {
			if (!tree->remove(vals[i])) {
				ERROR("remove failed!");
				return EXIT_FAILURE;
			}
			if (tree->find(vals[i]) != nullptr) {
				ERROR("Query failed after remove!");
				return EXIT_FAILURE;
			}
		}
	}

	PRINT("Test passed!!");

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
