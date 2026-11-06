# subopt-kaks top-level build
#
#   make                build squid (Sean Eddy's library) and the C tools
#   make test           run the regression tests against bundled data
#   make install        install C tools into $(PREFIX)/bin
#                       (PREFIX defaults to this directory, giving ./bin)
#   make install PREFIX=/usr/local
#   make rust           build the Rust port (rust/subopt-kaks, binary `kaks`)
#   make install-rust   cargo install kaks into $(PREFIX)/bin
#   make clean          remove build products (keeps squid's configure output)
#   make distclean      also remove squid's generated Makefile/config
#
# Any variable (CC, CFLAGS, PREFIX, DATADIR, DESTDIR) is passed to src/Makefile.

PREFIX ?= $(CURDIR)
SQUID   = squid
SRC     = src
RUST    = rust/subopt-kaks
CC     ?= gcc

export CC PREFIX

.PHONY: all squid tools test install uninstall rust test-rust install-rust \
        clean distclean

all: tools

# --- squid ---------------------------------------------------------------
$(SQUID)/Makefile: $(SQUID)/configure $(SQUID)/Makefile.in
	cd $(SQUID) && ./configure CC="$(CC)"

$(SQUID)/libsquid.a: $(SQUID)/Makefile
	$(MAKE) -C $(SQUID) libsquid.a

squid: $(SQUID)/libsquid.a

# --- C tools -------------------------------------------------------------
tools: $(SQUID)/libsquid.a
	$(MAKE) -C $(SRC)

test: tools
	$(MAKE) -C $(SRC) test

install: tools
	$(MAKE) -C $(SRC) install

uninstall:
	$(MAKE) -C $(SRC) uninstall

# --- Rust port -----------------------------------------------------------
rust:
	cd $(RUST) && cargo build --release

test-rust:
	cd $(RUST) && cargo test --release
	tst/run_tests.sh --rust $(CURDIR)/$(RUST)/target/release

install-rust:
	cargo install --path $(RUST) --root $(PREFIX)

# --- cleanup -------------------------------------------------------------
clean:
	$(MAKE) -C $(SRC) clean
	-[ -f $(SQUID)/Makefile ] && $(MAKE) -C $(SQUID) clean
	rm -f $(SQUID)/libsquid.a

distclean: clean
	-[ -f $(SQUID)/Makefile ] && $(MAKE) -C $(SQUID) distclean
	rm -rf bin
	-cd $(RUST) && cargo clean
