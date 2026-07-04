CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

../bin/sgclinalg: main.cpp sms.hpp gf2.hpp
	mkdir -p ../bin
	$(CXX) $(CXXFLAGS) -o $@ main.cpp

test: ../bin/sgclinalg
	bash test/run_tests.sh

clean:
	rm -f ../bin/sgclinalg

.PHONY: test clean
