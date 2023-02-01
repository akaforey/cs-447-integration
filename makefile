integrate: main.cpp
	g++ -mavx -pthread -O3 -Wall -o integrate main.cpp

debug: main.cpp
	g++ -mavx -pthread -g -Wall -o debug main.cpp

clean:
	rm integrate debug
