all: obj test
obj:
	g++ -std=c++14 -g -c pllml.cpp
test: obj
	g++ -std=c++14 -g -o test pllml.o test.cpp -lpll-sse3
clean:
	rm -f test *.o
