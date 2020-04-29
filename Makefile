PRG = Heat2D_MPI
EXECUTABLE = $(PRG).exe

#Use the GNU Compiler
MPI_CC = mpicc
#DEBUG = -g     %% uncomment this line for adding debugging option !Be careful: it slows down the execution process
CFLAGS = -O3 -std=c99 $(DEBUG)    # -O3 Optimization Level three for speed up
LIBS = -lm                        # -lm mathematical library
SOURCES = Therm2Dplate_mpi.c mpi_subroutines.c 
DEPS = parameters.h mpi_subroutines.h
OBJECTS = $(SOURCES:.c=.o)

%.o: %.c $(DEPS)
	 $(MPI_CC) -c -o $@ $< $(CFLAGS)

$(EXECUTABLE): $(OBJECTS)
	$(MPI_CC) $(LDFLAGS) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o $(EXECUTABLE) core
       
.PHONY: tar
tar:
	tar -cvf $(PRG).tar  $(SOURCES) $(DEPS) Makefile
	gzip $(PRG).tar

.PHONY: remove 

remove:
	rm -f processor*
