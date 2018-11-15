all:
	cd lib/gsl-2.4/ &&\
	mkdir installdir &&\
	./configure --prefix=${CURDIR}/lib/gsl-2.4/installdir &&\
	make && make install
	cd src && \
	gcc -Wall -I ../lib/gsl-2.4/installdir/include/ -c calc_deprat.c && \
	gcc calc_deprat.o ../lib/gsl-2.4/installdir/lib/libgsl.a ../lib/gsl-2.4/installdir/lib/libgslcblas.a -lm -o calc_deprat
	gcc -Wall -I ../lib/gsl-2.4/installdir/include/ -c calc_nclust.c && \
	gcc calc_nclust.o ../lib/gsl-2.4/installdir/lib/libgsl.a ../lib/gsl-2.4/installdir/lib/libgslcblas.a -lm -o calc_nclust
	mkdir ../bin && \
	cp calc_deprat calc_nclust ../bin
	cp script/extract_nclust.pl script/calc_cna.pl bin/
