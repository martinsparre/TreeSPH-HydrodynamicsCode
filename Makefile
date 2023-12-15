CC = gcc -pg -fopenmp
CFLAGS= -g -Wall  $(OPT) -O2
#"-pg" is for "gprof"
# gprof Tree gmon.out > profile_gcc_gravity.txt


#CC=icc -openmp
#CFLAGS= -g -O2 -marchcore2 -w1 $(OPT)
#-marchcore2: optimized for Intel Dual Core

#CFLAGS += -pedantic

LIBS= -lm

EXEC = TreeSPH
OBJS = all.o tree.o  main.o initialize.o print.o run.o density.o gravity.o sph.o
INCL = all.h Makefile

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) logfile.log
