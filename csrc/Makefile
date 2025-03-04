CC        = g++
NVCC      = nvcc
LD        = $(NVCC)
CFLAGS    = -Wall -O3 -I /usr/local/cuda/include
LIBS      = -lcufft -lcublas

#NVCCPARMS += -arch sm_13

PROG=cupoisson

OBJS =  main.o poisson.o utils.o

all: $(PROG)
       
$(PROG):	$(OBJS)
		$(LD) -o $@ $(OBJS) $(LIBS)
       

NVCCINC=-I $(CUDASDK)/common/inc

.SUFFIXES:

%.o:	%.c
		$(CC) $(CFLAGS) -o $@ -c $<

%.o:	%.cu
		$(NVCC) $(NVCCPARMS) -o $@ -c $<

clean:
	rm -f *.o $(PROG) *.linkinfo
