#! /bin/sh

INC1=/home/wyf/SeismicUnix/include
INC2=/home/wyf/SeismicUnix/src/Complex/include
LIK1=/home/wyf/SeismicUnix/lib
LIK2=/home/wyf/SeismicUnix/src/Complex/lib
LIB=-lsu -lpar -lcwp -lm

all: viscocpml.c
	g++ -o a.out viscocpml.c -fopenmp -I$(INC1) -I$(INC2) -L$(LIK1) -L$(LIK2)  $(LIB)
	clear        
	nohup ./a.out&
