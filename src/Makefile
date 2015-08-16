$cc=g++
all: act alt atm
act: act.o topiclib.o cokus.o util.o
	g++ -o act act.o topiclib.o cokus.o util.o -lm
alt: alt.o topiclib.o cokus.o util.o
	g++ -o alt alt.o topiclib.o cokus.o util.o -lm
atm: atm.o topiclib.o cokus.o util.o
	g++ -o atm atm.o topiclib.o cokus.o util.o -lm
act.o: act.cc
	g++ -c -g act.cc
alt.o: alt.cc
	g++ -c -g alt.cc
atm.o: atm.cc
	g++ -c -g atm.cc
topiclib.o: 
	g++ -c -g topiclib.cc
cokus.o:
	g++ -c -g cokus.cc
util.o:
	g++ -c -g util.cc
clean:
	rm -rf *.o
