include Makefile.inc

all: bottleneckBipartiteMatching_runner bvnGreedy mc64main doMatchingExps compareBisectionBasedOnMC64J3_Bottled_runner measureMMakerTime

#%.o: %.c
#	$(CC) -c $< $(CFLAGS)

bottleneckBipartiteMatching_runner: bottleneckBipartiteMatching_runner.c bottleneckBipartiteMatching.o mmio.o matrixUtils.o extern/cheap.o extern/matching.o
	make -C extern/
	$(CC)  $^ -o $@ $(CFLAGS) -lm

compareBisectionBasedOnMC64J3_Bottled_runner: compareBisectionBasedOnMC64J3_Bottled_runner.c bottleneckBipartiteMatching.o mmio.o matrixUtils.o  extern/cheap.o extern/matching.o
	make -C extern/
	$(CC)  $^ -o $@ $(CFLAGS) -lm

bvnGreedy: bvnGreedy.c bottleneckBipartiteMatching.o mmio.o matrixUtils.o extern/cheap.o extern/matching.o
	make -C extern/
	$(CC)  $^ -o $@ $(CFLAGS) -lm

mc64main: mc64main.c mmio.c matrixUtils.c mc64a.o
#	$(FCC) -c mc64a.F	 $(FFLAGS)
	$(CC)  $^ -o $@ $(CFLAGS) $(FLIBDIR) -lgfortran -lm

mc64a.o: mc64a.F
	$(FCC) -c mc64a.F	 $(FFLAGS)

doMatchingExps: doMatchingExps.c bottleneckBipartiteMatching.o mmio.o matrixUtils.o extern/cheap.o extern/matching.o
	$(FCC) -c mc64a.F	 $(FFLAGS)
	$(CC)  $^ -o $@ mc64a.o $(CFLAGS) $(FLIBDIR) -lgfortran -lm

measureMMakerTime: measureMMakerTime.c bottleneckBipartiteMatching.o mmio.o matrixUtils.o extern/cheap.o extern/matching.o
	$(CC)  $^ -o $@  $(CFLAGS) -lm

clean:
	make clean -C extern
	rm -f *.o bottleneckBipartiteMatching_runner compareBisectionBasedOnMC64J3_Bottled_runner bvnGreedy mc64main doMatchingExps measureMMakerTime

