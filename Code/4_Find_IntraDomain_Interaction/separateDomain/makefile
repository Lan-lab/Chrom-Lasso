OMPFLAG = -g -fopenmp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -s
CFLAGS = -I /srv/gsfs0/projects/pritchard/Xun/lib/include -I /home/xlan/Tools/GSL/gsl-1.16/ -I /usr/include/libxml2/ -Wall -c -g -fopenmp -std=c++11
LFLAGS = -L /srv/gsfs0/projects/pritchard/Xun/lib/lib -L /srv/gsfs0/projects/pritchard/Xun/lib/gsl-1.16/lib -L /srv/gsfs0/software/boost/1.51.0/lib/ -Wall $(OMPFLAG) -lmlpack -lboost_unit_test_framework -lboost_program_options -lm -std=c++11

separateDomain: separateDomain.o
	g++ $(LFLAGS) separateDomain.o -o separateDomain -Wl,-Bstatic -lboost_filesystem -lboost_system -Wl,-Bdynamic -g -fopenmp -l gsl -lgslcblas

separateDomain.o:
	gcc $(CFLAGS)  separateDomain.cpp -o separateDomain.o

clean:
	rm -rf *.o separateDomain
