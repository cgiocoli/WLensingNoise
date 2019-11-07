#####################################################
#                                                   #
#                        Makefile                   #
#                                                   #
#                         cgiocoli@gmail.com        #
#####################################################

# executable name
PROG = $(HOME)/bin/CreateWLNoise-0.9

PROG2 = $(HOME)/bin/CreateWLNoise_CoDECS_TEST_gsl

MAIN = noise_map.cpp

# .cpp lib internal in MOKA
SOURCES = noisekappa.cpp \
	  ../Moka/utilities.cpp \
	  ../WeakLMoka/power2D.cpp 


# gsl, cfitsio, CCfits, fftw
LIBS = -L/home/giocoli/lib/gsl-1.13/lib  -lgsl -lgslcblas \
       -L/home/giocoli/lib/cfitsio/lib \
       -L/home/giocoli/lib/CCfits/lib  -lCCfits -lcfitsio \
       -L/home/giocoli/lib/fftw-3.2.2/lib -lfftw3 -lm 

# gsl, cfitsio, CCfits, fftw  
ALLFLAGS = -I/home/giocoli/lib/gsl-1.13/include/gsl \
	   -I/home/giocoli/lib/gsl-1.13/include \
           -I/home/giocoli/lib/cfitsio/include \
           -I/home/giocoli/lib/CCfits/include \
           -I/home/giocoli/lib/fftw-3.2.2/include
# 
RELEASE = -O2 
# 
DEBUG = -g
# compiler  
#CC = g++
CC = /usr/bin/g++ -std=c++11 -Ofast -O3
#
RM = rm -f -r
#
OBJ = $(SOURCES:.cpp=.o)
#

CFLAGS=$(ALLFLAGS) $(DEBUG)
#
CLEAR = clear

default: ord
ord: 
	$(CC) -c $(SOURCES) $(CFLAGS)
	ar r libmokas.a *.o
	$(CC) $(MAIN) -L. -lmokas $(CFLAGS) -o $(PROG) $(LIBS) 

test: 
	$(CC) -c $(SOURCES) $(CFLAGS)
	ar r libmokas.a *.o
	$(CC) testgaussnoise.cpp -L. -lmokas $(CFLAGS) -o $(PROG2) $(LIBS) 

clean:
	$(RM) $(PROG) $(OBJ) $(PROG2)
