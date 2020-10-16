OMPFLAG = -g -fopenmp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -s
CFLAGS = -I /srv/gsfs0/projects/pritchard/Xun/lib/include -I /home/xlan/Tools/GSL/gsl-1.16/ -I ./ -I ../ -I /usr/include/libxml2/ -Wall -c -g -fopenmp -std=c++11
LFLAGS = -L /srv/gsfs0/projects/pritchard/Xun/lib/lib -L /srv/gsfs0/projects/pritchard/Xun/lib/gsl-1.16/lib -L /srv/gsfs0/software/boost/1.51.0/lib/ -Wall $(OMPFLAG) -lmlpack -lboost_unit_test_framework -lboost_program_options -lm -std=c++11

findIntraDomainInteraction: findIntraDomainInteraction.o IRLS.o LogLink.o
	g++ $(LFLAGS) findIntraDomainInteraction.o IRLS.o LogLink.o -o findIntraDomainInteraction -Wl,-Bstatic -lboost_filesystem -lboost_system -Wl,-Bdynamic -g -fopenmp -l gsl -lgslcblas

findIntraDomainInteraction.o: ../IRLS_glm/IRLS.h
	gcc $(CFLAGS)  findIntraDomainInteraction.cpp -o findIntraDomainInteraction.o

IRLS.o: ../IRLS_glm/IRLS.h
	gcc $(CFLAGS) ../IRLS_glm/IRLS.cc -o IRLS.o

LogLink.o: ../IRLS_glm/LogLink.h ../IRLS_glm/LinkFunc.h
	gcc $(CFLAGS) ../IRLS_glm/LogLink.cc -o LogLink.o


clean:
	rm -rf *.o findIntraDomainInteraction
