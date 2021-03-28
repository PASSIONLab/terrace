#!/bin/python

import sys

infile = sys.argv[1]
out_queries = sys.argv[2]
out_inserts = sys.argv[3]
outparts = infile.split("_")
threads = outparts[-2]
src_vtx = outparts[-1]
graph = outparts[0]
rounds= 10 # if you need to change this later
with open(infile, 'r') as inf:
  with open(out_queries, 'a') as outq_f:
    with open(out_inserts, 'a') as outi_f:
      for line in inf:
        # deal with queries
        if line.startswith("AVG"):
          # parse it
          parts = line.split()
          time = parts[2].split("=")[-1]
          alg = parts[1].split("=")[-1]
          # structure,algorithm,src,graph,threads,time
          outq_f.write("gc," + str(rounds) + "," + alg + "," + src_vtx + "," + graph + "," + threads + "," + time + "\n")

        # parse batch inserts
        if line.startswith("batch_size"):
          parts = line.split(",")
          update_type = parts[1].split()[1][:-1] 
          batch_size = parts[0].split()[-1]
          throughput = parts[2].split()[-1]
          # structure,update type,batch size,graph,throughput
          outi_f.write("gc," + update_type + "," + batch_size + "," + graph + "," + throughput + "\n")
