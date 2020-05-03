INSTALLDIR=/usr/local/bin

src = $(wildcard *.cpp *.c *.cxx)
csrc = $(wildcard *.c)
obj_cpp = $(src:.cpp=.o)
obj = $(obj_cpp:.c=.o)
hdr = $(wildcard *.h *.hpp)

CC=gcc
COPT=-g -lm
COPT2=-g -lm

#COPT=-O3 -lm
#COPT2=-O3 -lm -DGERBER_TEST

CXX=g++
CXXOPT = -g --std=c++11
#CXXOPT = -O3 --std=c++11
CXXLIB = -lstdc++ -lboost_thread -lm

CCFN = $(wildcard *.c)
CPPFN = $(wildcard *.cpp)

%.o: %.c
	$(CC) $(COPT) -c $< -o $@ -lm

%.o: %.cpp
	$(CXX) $(CXXOPT) -c $< -o $@

%.o: %.cxx
	$(CC) $(COPT) -c $< -o $@ -lm

survey_prop: $(obj) $(hdr)
	$(CXX) $(CXXOPT) -o $@ $(obj) $(CXXLIBS)

test_survey_prop: $(csrc)
	$(CC) $(COPT2) -o $@ $^ $(COPT2)

.PHONY: debug
debug: $(obj) $(hdr)
	$(CXX) -g -o survey_prop_debug $^ $(CXXLIBS)

.PHONY: clean
clean:
	rm -f *.o survey_prop
	rm -f *.o test_survey_prop

.PHONY: install
install: survey_prop
	mkdir -p $(INSTALLDIR)
	cp $< $(INSTALLDIR)
