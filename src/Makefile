# subopt-kaks: suboptimal alignments and Ka/Ks estimation
#
# Normal use is through the top-level Makefile (which also builds squid):
#     make                 build everything
#     make install         install into ./bin (PREFIX defaults to the repo)
#     make install PREFIX=/usr/local
#
# Variables you may override on the command line:
#     CC, CFLAGS, LDFLAGS  compiler and flags
#     PREFIX               install root (default: repository root)
#     DATADIR              where matricies/ lives at run time
#                          (default: $(PREFIX)/share/subopt-kaks, or the
#                          repository root when PREFIX is the repository)
#     DESTDIR              staging prefix for packaging
#
# The run-time matrix directory is compiled into util.o as SUBOPT_DATADIR.
# A `.datadir` stamp forces util.o to rebuild whenever DATADIR changes, so
# `make install PREFIX=...` always installs binaries that know their data.

ROOT      := $(abspath ..)
PREFIX    ?= $(ROOT)
ifeq ($(abspath $(PREFIX)),$(ROOT))
DATADIR   ?= $(ROOT)
MANDIR    ?= $(ROOT)/man
else
DATADIR   ?= $(PREFIX)/share/subopt-kaks
MANDIR    ?= $(PREFIX)/share/man/man1
endif
BINDIR    ?= $(PREFIX)/bin
DESTDIR   ?=

CC        ?= gcc
SQUIDDIR  ?= $(ROOT)/squid
SQUIDLIB  ?= $(SQUIDDIR)/libsquid.a
CFLAGS    ?= -O2 -g
CPPFLAGS  += -I$(SQUIDDIR) -DSUBOPT_DATADIR='"$(DATADIR)"'
LDFLAGS   ?=
LIBS       = $(SQUIDLIB) -lm

INSTALL      ?= install
INSTALL_PROG ?= $(INSTALL) -m 755
INSTALL_DATA ?= $(INSTALL) -m 644

PROGS     = zuker_align cds_zuker_align yn00_zuker_align yn00_cds_optimal \
            yn00_cds_prealigned
MANPAGES  = $(addprefix $(ROOT)/man/,$(addsuffix .1,$(PROGS)))

ALIGN_OBJS = zuker.o util.o
CDS_OBJS   = $(ALIGN_OBJS) mrtrans.o
YN00_OBJS  = yn00_simple.o tools.o
PROG_OBJS  = $(addsuffix .o,$(PROGS))

.PHONY: all install uninstall clean test

all: $(PROGS)

zuker_align:         $(ALIGN_OBJS) zuker_align.o
cds_zuker_align:     $(CDS_OBJS) cds_zuker_align.o
yn00_zuker_align:    $(CDS_OBJS) $(YN00_OBJS) yn00_zuker_align.o
yn00_cds_optimal:    $(CDS_OBJS) $(YN00_OBJS) yn00_cds_optimal.o
yn00_cds_prealigned: util.o $(YN00_OBJS) yn00_cds_prealigned.o

$(PROGS):
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

# header dependencies
zuker.o util.o zuker_align.o: align.h util.h
mrtrans.o cds_zuker_align.o: mrtrans.h util.h
yn00_simple.o tools.o: paml.h treesub.c util.h
yn00_zuker_align.o yn00_cds_optimal.o: align.h mrtrans.h paml.h util.h
yn00_cds_prealigned.o: align.h mrtrans.h util.h

# rebuild util.o when the compiled-in data directory changes
.datadir: FORCE
	@echo '$(DATADIR)' | cmp -s - $@ || echo '$(DATADIR)' > $@
util.o: .datadir
FORCE:

$(SQUIDLIB):
	@echo "squid library $(SQUIDLIB) not found; run 'make' at the top level" >&2
	@exit 1

install: $(PROGS)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL_PROG) $(PROGS) $(DESTDIR)$(BINDIR)
ifneq ($(abspath $(MANDIR)),$(ROOT)/man)
	$(INSTALL) -d $(DESTDIR)$(MANDIR)
	$(INSTALL_DATA) $(MANPAGES) $(DESTDIR)$(MANDIR)
endif
ifneq ($(abspath $(DATADIR)),$(ROOT))
	$(INSTALL) -d $(DESTDIR)$(DATADIR)/matricies/aa $(DESTDIR)$(DATADIR)/matricies/nt
	$(INSTALL_DATA) $(ROOT)/matricies/README $(DESTDIR)$(DATADIR)/matricies
	$(INSTALL_DATA) $(ROOT)/matricies/aa/* $(DESTDIR)$(DATADIR)/matricies/aa
	$(INSTALL_DATA) $(ROOT)/matricies/nt/* $(DESTDIR)$(DATADIR)/matricies/nt
endif

uninstall:
	rm -f $(addprefix $(DESTDIR)$(BINDIR)/,$(PROGS))
ifneq ($(abspath $(MANDIR)),$(ROOT)/man)
	rm -f $(addprefix $(DESTDIR)$(MANDIR)/,$(notdir $(MANPAGES)))
endif
ifneq ($(abspath $(DATADIR)),$(ROOT))
	rm -rf $(DESTDIR)$(DATADIR)/matricies
endif

test: $(PROGS)
	$(ROOT)/tst/run_tests.sh $(CURDIR)

clean:
	rm -f *~ *.o $(PROGS) .datadir
