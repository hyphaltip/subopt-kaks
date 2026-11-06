# kaks (Rust port of subopt-kaks)
#
#   make                 cargo build --release  -> target/release/kaks
#   make test            cargo test + regression tests against ../../tst
#   make install         copy the built binary into $(BINDIR)
#   make install PREFIX=/usr/local      (-> /usr/local/bin/kaks)
#   make install BINDIR=$HOME/bin
#   make uninstall
#   make clean
#
# DESTDIR is honoured for staged installs. The binary is fully static with
# respect to this package (BLOSUM62 and the genetic codes are embedded), so
# nothing else needs to be installed.

PREFIX  ?= $(HOME)/.local
BINDIR  ?= $(PREFIX)/bin
DESTDIR ?=
CARGO   ?= cargo
INSTALL ?= install

BIN     = kaks
TARGET  = target/release/$(BIN)
TST     = $(abspath ../../tst)

.PHONY: all build test install uninstall clean

all: build

build $(TARGET):
	$(CARGO) build --release

test: build
	$(CARGO) test --release
	$(TST)/run_tests.sh --rust $(abspath target/release)

install: build
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 755 $(TARGET) $(DESTDIR)$(BINDIR)/$(BIN)
	@echo "installed $(DESTDIR)$(BINDIR)/$(BIN)"

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(BIN)

clean:
	$(CARGO) clean
