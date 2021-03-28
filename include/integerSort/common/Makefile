include parallelDefs

CHECK = isortCheck
COMMON = IO.h sequenceIO.h parseCommandLine.h parallel.h blockRadixSort.h transpose.h runTests.py sequence.h utils.h gettime.h


$(CHECK): $(CHECK).C $(COMMON)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $(CHECK) $(CHECK).C

$(COMMON) :
	ln -s ../../common/$@ .

clean :
	rm -f $(CHECK)

cleansrc :
	make -s clean
	rm -f $(COMMON) *.pyc
