# kaks (Rust port of subopt-kaks)
#
#   make                 cargo build --release  -> target/release/kaks
#   make test            cargo test + regression tests against ../../tst
#   make install         copy the built binary into $(BINDIR)
#   make install PREFIX=/usr/local      (-> /usr/local/bin/kaks)
#   make install BINDIR=$HOME/bin
#   make uninstall
#   make wasi            WASI command-line build -> target/wasm32-wasip1/release/kaks.wasm
#   make web             browser build -> ../../web/pkg (needs wasm-bindgen-cli 0.2.130)
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

.PHONY: all build test install uninstall clean wasi web

WEB     = $(abspath ../../web)
WASM    = target/wasm32-unknown-unknown/release/subopt_kaks.wasm

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

# Runs anywhere with a WASI runtime (wasmtime, wasmer, node --experimental-wasi).
wasi:
	rustup target add wasm32-wasip1
	$(CARGO) build --release --target wasm32-wasip1
	@echo "built target/wasm32-wasip1/release/kaks.wasm"

# Browser bundle: wasm + JS glue into $(WEB)/pkg, next to $(WEB)/index.html.
# wasm-bindgen-cli must match the wasm-bindgen crate version in Cargo.lock:
#     cargo install wasm-bindgen-cli --version 0.2.130
web:
	rustup target add wasm32-unknown-unknown
	$(CARGO) build --release --target wasm32-unknown-unknown --no-default-features --features wasm
	wasm-bindgen --target web --no-typescript --out-dir $(WEB)/pkg $(WASM)
	@echo "built $(WEB)/pkg; serve $(WEB) over HTTP (make -C ../.. serve-web)"

clean:
	$(CARGO) clean
