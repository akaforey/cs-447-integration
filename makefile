run: integrate
	./integrate 0 1 100 10

integrate: main.cpp
	g++ -pthread -g -Wall -o integrate main.cpp

clean:
	rm integrate
