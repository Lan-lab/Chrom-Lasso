OMPFLAG = -g -fopenmp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -s
#CFLAGS = -I /Share/home/lanxun4/yaoli/projects/lujz/include -I /home/xlan/Tools/GSL/gsl-1.16/ -I ./ -I ../ -I /usr/include/libxml2/ -Wall -c -g -fopenmp -std=c++11
#LFLAGS = -L /Share/home/lanxun4/yaoli/projects/lujz/lib -L /Share/home/lanxun4/yaoli/projects/lujz/lib64 -L /srv/gsfs0/projects/pritchard/Xun/lib/gsl-1.16/lib -L /Share/home/lanxun4/yaoli/toolshed/boost_1_51_0/stage/lib  -L /srv/gsfs0/projects/pritchard/Xun/lib/gsl-1.16/lib -Wall $(OMPFLAG) -lmlpack -lboost_unit_test_framework -lboost_program_options -lm -std=c++11
CFLAGS = -I ./ -I ../ -I /usr/include/libxml2 -I /Share/home/YaoLi/E_only/include/scythestat -Wall -c -g -fopenmp -std=c++11
LFLAGS = -L /Share/home/YaoLi/E_only/dependencies/boost_1_51_0/stage/lib -L /Share/home/YaoLi/E_only/lib64 -L /Share/home/YaoLi/E_only/lib -Wall -lmlpack -lboost_unit_test_framework -lboost_program_options -lm -std=c++11
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
