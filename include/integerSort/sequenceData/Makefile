include parallelDefs

COMMON = dataGen.h utils.h sequenceIO.h IO.h parallel.h sequence.h parseCommandLine.h
SEQUENCEGEN = $(COMMON) 
GENERATORS = randomSeq almostSortedSeq exptSeq trigramSeq addDataSeq trigramString

.PHONY: all clean
all: $(GENERATORS)

$(COMMON) :
	ln -s ../../common/$@ .

addDataSeq : addDataSeq.C sequenceData.h $(SEQUENCEGEN)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ $@.C

randomSeq : randomSeq.C sequenceData.h $(SEQUENCEGEN)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ $@.C

almostSortedSeq : almostSortedSeq.C sequenceData.h $(SEQUENCEGEN)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ $@.C

exptSeq : exptSeq.C sequenceData.h $(SEQUENCEGEN)
	$(PCC) $(PCFLAGS) $(PLFLAGS) -o $@ $@.C

trigrams.o : trigrams.C $(SEQUENCEGEN) sequence.h
	$(PCC) $(PCFLAGS) -c trigrams.C

trigramSeq.o : trigramSeq.C $(SEQUENCEGEN)
	$(PCC) $(PCFLAGS) -c trigramSeq.C

trigramSeq : trigramSeq.o trigrams.o 
	$(PCC) $(PLFLAGS) -o $@ $@.o trigrams.o

trigramString.o : trigramString.C $(SEQUENCEGEN)
	$(PCC) $(PCFLAGS) -c trigramString.C

trigramString : trigramString.o trigrams.o 
	$(PCC) $(PLFLAGS) -o $@ $@.o trigrams.o

clean :
	rm -f *.o $(GENERATORS)
	make clean -s -C data

cleansrc : 
	make -s clean
	rm -f $(COMMON) 
