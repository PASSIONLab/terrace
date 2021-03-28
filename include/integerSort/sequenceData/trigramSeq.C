#include "parallel.h"
#include "sequenceIO.h"
#include "parseCommandLine.h"
using namespace benchIO;

char** trigramWords(intT s, intT e);

int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"<size> <outfile>");
  pair<intT,char*> in = P.sizeAndFileName();
  char **A = trigramWords(0,in.first);
  return writeSequenceToFile(A,in.first,in.second);
}
