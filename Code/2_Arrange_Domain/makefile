
OMPFLAG = -g -fopenmp -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -s
CFLAGS =  -I /Share/app/R-20180408/R3.4.4/lib64/R/include  -I /Share/app/R-20180408/R3.4.4/lib64/R/library/Rcpp/include -I /Share/home/lanxun4/miniconda3/include/ -I /Share/home/lanxun/R/x86_64-pc-linux-gnu-library/3.4/RInside/include -I /Share/home/lanxun4/yaoli/projects/lujz/include  -Wall $(CXXFLAGS) -c -g -fopenmp
# LFLAGS = -L /Share/home/YaoLi/E_only/lib64 -L /Share/home/YaoLi/E_only/lib -L /Share/home/lanxun4/yaoli/projects/lujz/lib -L /Share/home/lanxun4/yaoli/projects/lujz/lib64 -L Share/home/lanxun/Xun/scg48U/home/Tools/GSL/gsl-1.16/lib -L /Share/home/lanxun/Xun/scg4BU/Tools/R-3.0.1/lib -L /Share/home/lanxun/Xun/scg4BU/Tools/R-3.0.1/library/Rcpp/lib -L /Share/home/lanxun/Xun/scg4BU/Tools/R-3.0.1/library/RInside/lib -L /Share/home/YaoLi/E_only/dependencies/boost_1_51_0/stage/lib -L /Share/home/YaoLi/E_only/lib64 -L /Share/home/YaoLi/E_only/lib -Wall $(OMPFLAG) $(LDLIBS) -lm -lR -lRInside -lreadline -lgsl -lRcpp -lRblas -lRlapack -L/Share/app/gcc/4.9.4/lib64 -lgfortran
LFLAGS = -L /Share/home/YaoLi/E_only/lib64 -L /Share/home/YaoLi/E_only/lib -L /Share/home/lanxun4/yaoli/projects/lujz/lib -L /Share/home/lanxun4/yaoli/projects/lujz/lib64 -L Share/home/lanxun/Xun/scg48U/home/Tools/GSL/gsl-1.16/lib -L /Share/app/R-20180408/R3.4.4/lib64/R/lib -L /Share/app/R-20180408/R3.4.4/lib64/R/library/Rcpp/lib -L /Share/home/lanxun/R/x86_64-pc-linux-gnu-library/3.4/RInside/lib -L /Share/home/YaoLi/E_only/dependencies/boost_1_51_0/stage/lib -L /Share/home/YaoLi/E_only/lib64 -L /Share/home/YaoLi/E_only/lib -Wall $(OMPFLAG) $(LDLIBS) -lm -lR -lRInside -lreadline -lgsl -lRblas -lRlapack -lgfortran

HiC_mixturePLD_singleThread: HiC_mixturePLD.o HiC_BOSS.o HiC_POSTDOC.o HiC_regionDetail.o IRLS.o LogLink.o
	g++ $(LFLAGS) HiC_mixturePLD.o HiC_BOSS.o HiC_POSTDOC.o HiC_regionDetail.o IRLS.o LogLink.o -o HiC_mixturePLD_singleThread -Wl,-Bstatic -lboost_filesystem -lboost_system -Wl,-Bdynamic -g -fopenmp -l gsl -lgslcblas

HiC_mixturePLD.o: HiC_mixturePLD_boss.h HiC_mixturePLD_postdoc.h
	g++ $(CFLAGS) HiC_mixturePLD_main.cpp -o HiC_mixturePLD.o

HiC_BOSS.o: HiC_mixturePLD_boss.h bothEndsMappedFragInfo_withCuttingSite.hpp regionDetail2D_pairSite.hpp
	g++ $(CFLAGS)  HiC_mixturePLD_boss.cpp -o HiC_BOSS.o

HiC_POSTDOC.o: HiC_mixturePLD_postdoc.h HiC_mixturePLD_boss.h
	g++ $(CFLAGS)  HiC_mixturePLD_postdoc.cpp -o HiC_POSTDOC.o

HiC_regionDetail.o: regionDetail2D_pairSite.hpp
	g++ $(CFLAGS) regionDetail2D_pairSite.cpp -o HiC_regionDetail.o

IRLS.o: ../IRLS_glm/IRLS.h
	g++ $(CFLAGS) ../IRLS_glm/IRLS.cc -o IRLS.o

LogLink.o: ../IRLS_glm/LogLink.h ../IRLS_glm/LinkFunc.h
	g++ $(CFLAGS) ../IRLS_glm/LogLink.cc -o LogLink.o

clean:
	rm -rf *.o HiC_mixturePLD

