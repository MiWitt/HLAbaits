

CXX      = g++   # for Linux RedHat 6.1, g++ version 2.95.2
CXXFLAGS = -O2 -Wall -I.


default: HLAbaits

clean:
	@echo "Cleaning up."
	@rm -f *.o HLAbaits

CMyException.o: CMyException.cpp CMyException.h
	${CXX} ${CXXFLAGS} -c -o CMyException.o CMyException.cpp

CMultiFasta.o: CMultiFasta.cpp CMultiFasta.h
	${CXX} ${CXXFLAGS} -c -o CMultiFasta.o CMultiFasta.cpp

main.o: main.cpp CMyException.o CMultiFasta.o
	${CXX} ${CXXFLAGS} -w -c -o main.o main.cpp
	
HLAbaits: main.o
	${CXX} -o HLAbaits main.o CMultiFasta.o CMyException.o
	



