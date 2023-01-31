integrate: main.cpp
	g++ -pthread -g -Wall -o integrate main.cpp

clean:
	rm integrate
