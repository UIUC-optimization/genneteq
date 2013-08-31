CC = g++
CFLAGS = -O3 -I/usr/local/include
DFLAGS = -g -pg -O3
LDFLAGS = -lm -L/usr/local/lib -lgsl -lgslcblas 

PROG = genneteq
DEBUG_PROG = genneteq-debug
HDRS = main.h matrix.h graph.h variables.h trees.h genneteq.h
SRCS = main.cpp graph.cpp trees.cpp genneteq.cpp gneInit.cpp gneFlow.cpp gnePotential.cpp gnePivot.cpp gnePivotRules.cpp #gneDebug.cpp

## This incantation says that the object files
## have the same name as the .c files, but with .o
OBJS = $(SRCS:.cpp=.o)
DOBJS = $(SRCS:.cpp=.d)

## This is the first rule (the default)
default : $(PROG)
debug : $(DEBUG_PROG)
all : $(PROG) $(DEBUG_PROG)

## Build the program from the .o's
$(PROG) : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROG) $(LDFLAGS)

$(DEBUG_PROG) : $(DOBJS) 
	$(CC) $(DFLAGS) $(DOBJS) -o $(DEBUG_PROG) $(LDFLAGS)

## Rules for the source files
.cpp.o: 
	$(CC) $(CFLAGS) -c -o $@ $<

%.d : %.cpp
	$(CC) $(DFLAGS) -c -o $@ $<

.PHONY: clean install

## Remove all the compilation and debugging files
clean :
	rm -f core $(OBJS) $(DOBJS)

install: $(PROG)
	cp $(PROG) ../bin/.

