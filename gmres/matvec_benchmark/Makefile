FLAGS = -Wall -Wextra
CC = gcc $(FLAGS)

all: matvecBenchmark.cpp
	g++-4.6 matvecBenchmark.cpp -o matvecBenchmark -L/opt/acml5.3.1/gfortran64/lib /opt/acml5.3.1/gfortran64/lib/libacml.a -lm -ldl -lpthread -lgfortran -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -g -O2
