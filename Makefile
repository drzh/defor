all:
	cd lib/gsl-2.4/ &&\
	mkdir installdir &&\
	./configure --prefix=${CURDIR}/installdir &&\
	make install
	cd src && \
	gcc -Wall -I ../lib/gsl-2.4/installdir/include/ -c calc_deprat.c && \
	gcc -L ../lib/gsl-2.4/installdir/lib/ -static calc_deprat.o -lgsl -lgslcblas -lm -o calc_deprat && \
	gcc -Wall -I ../lib/gsl-2.4/installdir/include/ -c calc_nclust.c && \
	gcc -L ../lib/gsl-2.4/installdir/lib/ -static calc_nclust.o -lgsl -lgslcblas -lm -o calc_nclust && \
	mkdir ../bin && \
	cp calc_deprat calc_nclust ../bin
