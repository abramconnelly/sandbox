FLAGS=-Wall -Winline -O3 -DQUEUE
#FLAGS=-Wall -Winline -O3 -funroll-loops -DFAST_ITERATION
#FLAGS=-Wall -Winline -O3 -funroll-loops -DQUEUE
#FLAGS=-g -Wall -pg
#FLAGS=-g
#CC=gcc

all:	sp merge verify walksat
sp: sp.c sp.h formula.h formula.o random.o queue.o
	${CC} ${FLAGS} sp.c random.o formula.o queue.o -o sp -lm
formula.o: formula.h formula.c
	${CC} ${FLAGS} -c formula.c -o formula.o
random.o: random.c random.h
	${CC} ${FLAGS} -c random.c -o random.o
queue.o: queue.c queue.h
	${CC} ${FLAGS} -c queue.c -o queue.o
clean:
	rm -f core *~ *.o  gmon.out *.tmp.* *.out *.dat
binclean: clean
	rm -f walksat sp spprof verify merge
walksat: walksat.c
	${CC} -lm -O6 walksat.c -o walksat
merge: merge.c
	${CC} -Wall merge.c -o merge
verify: verify.c
	${CC} -Wall verify.c -o verify
