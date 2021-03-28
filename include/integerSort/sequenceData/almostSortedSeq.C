#include "sequenceData.h"
#include "sequenceIO.h"
#include "parseCommandLine.h"
using namespace dataGen;
using namespace benchIO;

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-r <swaps>] [-t {int,double}] <size> <outfile>");
  pair<intT,char*> in = P.sizeAndFileName();
  elementType tp = elementTypeFromString(P.getOptionValue("-t","int"));
  intT n = in.first;
  char* fname = in.second;
  intT swaps = P.getOptionIntValue("-r",floor(sqrt((float) n)));
  
  switch(tp) {
  case intType: 
    return writeSequenceToFile(almostSorted<intT>(0, n, swaps), n, fname);
  case doubleT: 
    return writeSequenceToFile(almostSorted<double>(0, n, swaps), n, fname);
  default:
    cout << "genSeqRand: not a valid type" << endl;
  }
}
