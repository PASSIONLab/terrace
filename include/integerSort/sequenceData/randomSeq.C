#include "sequenceData.h"
#include "sequenceIO.h"
#include "parseCommandLine.h"
using namespace dataGen;
using namespace benchIO;

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"[-r <range>] [-t {int,double}] <size> <outfile>");
  pair<intT,char*> in = P.sizeAndFileName();
  elementType tp = elementTypeFromString(P.getOptionValue("-t","int"));
  intT n = in.first;
  char* fname = in.second;
  intT r = P.getOptionLongValue("-r",n);

  switch(tp) {
  case intType: return writeSequenceToFile(randIntRange(0,n,r), n, fname);
  case doubleT: return writeSequenceToFile(rand<double>(0, n), n, fname);
  default: cout << "genSeqRand: not a valid type" << endl;
  }
}
