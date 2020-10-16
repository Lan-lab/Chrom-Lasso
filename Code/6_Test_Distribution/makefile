OMPFLAG = -g -fopenmp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -s
CFLAGS = -I /Share/home/YaoLi/E_only/include -I ./ -I ../ -I /usr/include/libxml2/ -Wall -c -g -fopenmp -std=c++11
LFLAGS = -L /Share/home/YaoLi/E_only/lib -L /Share/home/YaoLi/E_only/lib64 -L /Share/home/YaoLi/E_only/dependencies/boost_1_51_0/stage/lib -Wall $(OMPFLAG) -lmlpack -lboost_unit_test_framework -lboost_program_options -lm -std=c++11

outputTestPosPvalue: outputTestPosPvalue.o IRLS.o LogLink.o
	g++ $(LFLAGS) outputTestPosPvalue.o IRLS.o LogLink.o -o outputTestPosPvalue -Wl,-Bstatic -lboost_filesystem -lboost_system -Wl,-Bdynamic -g -fopenmp -l gsl -lgslcblas

outputTestPosPvalue.o: ../IRLS_glm/IRLS.h
	gcc $(CFLAGS)  outputTestPosPvalue.cpp -o outputTestPosPvalue.o

IRLS.o: ../IRLS_glm/IRLS.h
	gcc $(CFLAGS) ../IRLS_glm/IRLS.cc -o IRLS.o

LogLink.o: ../IRLS_glm/LogLink.h ../IRLS_glm/LinkFunc.h
	gcc $(CFLAGS) ../IRLS_glm/LogLink.cc -o LogLink.o


clean:
	rm -rf *.o outputTestPosPvalue
