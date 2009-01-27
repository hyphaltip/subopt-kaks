CC        = gcc
# Intel compiler
#CC        = icc
COPY      = cp
MKDIR	  = mkdir -p
INSTALLDIR=../bin

PROGS     = zuker_align cds_zuker_align yn00_zuker_align yn00_cds_optimal \
	    yn00_cds_prealigned
OBJS      = zuker.o util.o
CDSOBJS   = $(OBJS) mrtrans.o 
YN00OBJS  = $(CDSOBJS) yn00_simple.o tools.o
CDS_OPT   = $(YN00OBJS)
PROG_OBJS = zuker_align.o cds_zuker_align.o yn00_zuker_align.o \
	    yn00_cds_optimal.o

SQUIDDIR  = ../squid
SQUIDINC  = $(SQUIDDIR)
SQUIDLIB  = $(SQUIDDIR)

#CFLAGS    = -g -DMALLOC -DDMALLOC_FUNC_CHECK -I$(SQUIDINC)
#CFLAGS    = -g  -I$(SQUIDINC) -O4
#LIBS      = -L$(SQUIDLIB) -lsquid -lm -ldmalloc
# Turn compiler optimization flags on
# athlon options
#CFLAGS    = -I$(SQUIDINC) -O4 -funroll-loops -fomit-frame-pointer \
#	    -march=athlon -mcpu=athlon 
# generic OSX options
CFLAGS = -I$(SQUIDINC) -O4 -funroll-loops -fomit-frame-pointer
# Intel options
#CFLAGS    = -I$(SQUIDINC) -O3 -tpp7 -unroll -axKW
LIBS      = -L$(SQUIDLIB) -lsquid -lm


all: $(PROGS)

install: $(PROGS)
	$(MKDIR) $(INSTALLDIR)
	$(COPY) $(PROGS) $(INSTALLDIR)

yn00_zuker_align: $(YN00OBJS) yn00_zuker_align.o
	$(CC) $(CFLAGS) -o $@ $(YN00OBJS) yn00_zuker_align.o $(LIBS)

zuker_align: $(OBJS) zuker_align.o
	$(CC) $(CFLAGS) -o $@ $(OBJS) zuker_align.o $(LIBS)

cds_zuker_align: $(CDSOBJS) cds_zuker_align.o
	$(CC) $(CFLAGS) -o $@ $(CDSOBJS) cds_zuker_align.o $(LIBS)

yn00_cds_optimal: $(CDS_OPT) yn00_cds_optimal.o
	$(CC) $(CFLAGS) -o $@ $(CDS_OPT) yn00_cds_optimal.o $(LIBS)

yn00_cds_prealigned: $(YN00OBJS) yn00_cds_prealigned.o
	$(CC) $(CFLAGS) -o $@ tools.o util.o yn00_simple.o yn00_cds_prealigned.o $(LIBS)


.c.o:
	$(CC) $(CFLAGS) -c $<
.c:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *~ $(OBJS) $(CDSOBJS) $(YN00OBJS) $(PROGS) $(PROG_OBJS)

zuker.o: align.h util.h

util.o: util.h

mrtrans.o: mrtrans.h mrtrans.c

zuker_align.o: zuker.o zuker_align.c

cds_zuker_align.o: zuker.o mrtrans.o cds_zuker_align.c

suboptaligner: suboptaligner.o align.o 
	$(CC) $(CFLAGS) -DDEBUG -o $@ suboptaligner.o align.o $(LIBS)

