.DEFAULT_GOAL := all

CC ?= cc
AR ?= ar
RANLIB ?= ranlib
MAKE ?= make
BUILD_DIR ?= build

HTSLIB_VERSION ?= 1.23.1
LMDB_VERSION ?= LMDB_0.9.31

DOWNLOAD_DIR := $(BUILD_DIR)/downloads
VENDOR_SRC_DIR := $(BUILD_DIR)/vendor/src
OBJ_DIR := $(BUILD_DIR)/obj

HTSLIB_URL ?= https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2
LMDB_URL ?= https://github.com/LMDB/lmdb/archive/$(LMDB_VERSION).tar.gz

HTSLIB_ARCHIVE := $(DOWNLOAD_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
LMDB_ARCHIVE := $(DOWNLOAD_DIR)/lmdb-$(LMDB_VERSION).tar.gz

HTSLIB_SRC := $(VENDOR_SRC_DIR)/htslib-$(HTSLIB_VERSION)
LMDB_SRC := $(VENDOR_SRC_DIR)/lmdb-$(LMDB_VERSION)/libraries/liblmdb

HTSLIB_PREFIX := $(abspath $(BUILD_DIR)/vendor/htslib)
LMDB_PREFIX := $(abspath $(BUILD_DIR)/vendor/lmdb)

HTSLIB_LIB := $(HTSLIB_PREFIX)/lib/libhts.a
LMDB_LIB := $(LMDB_PREFIX)/lib/liblmdb.a

HTSLIB_CONFIGURE_FLAGS ?= --disable-libcurl --disable-plugins --disable-s3 --disable-gcs --without-libdeflate
SYS_LIBS ?= -lz -lbz2 -llzma -lm -lpthread
EXE_LDFLAGS ?= -no-pie

WARNFLAGS := -Wall -Wextra -Wno-unused-parameter -Wshadow
CFLAGS ?= -O3 -g
CPPFLAGS += -Iinclude -Isrc -I$(HTSLIB_PREFIX)/include -I$(LMDB_PREFIX)/include \
	-D_POSIX_C_SOURCE=200809L -D_GNU_SOURCE
LDFLAGS +=

ISAM_SOURCES := \
	src/common.c \
	src/backend_lmdb.c \
	src/hash.c \
	src/rebuild.c \
	src/store.c \
	src/view.c \
	src/cmd_index.c \
	src/cmd_import.c \
	src/cmd_sort.c \
	src/cmd_view.c \
	src/cmd_annotate.c \
	src/cmd_update.c \
	src/cmd_rehint.c \
	src/api.c \
	src/htsisam_dispatch.c

ISAM_OBJECTS := $(patsubst src/%.c,$(OBJ_DIR)/src/%.o,$(ISAM_SOURCES))
MAIN_OBJECT := $(OBJ_DIR)/src/main.o

all: $(BUILD_DIR)/libisam.a $(BUILD_DIR)/isamtools

deps: $(HTSLIB_LIB) $(LMDB_LIB)

test: all
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/functional.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/lmdb_store.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/lmdb_alignment.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/lmdb_compat.sh

clean:
	rm -rf $(OBJ_DIR) $(BUILD_DIR)/libisam.a $(BUILD_DIR)/isamtools

distclean:
	rm -rf $(BUILD_DIR)

print-config:
	@printf 'BUILD_DIR=%s\n' '$(BUILD_DIR)'
	@printf 'HTSLIB_VERSION=%s\n' '$(HTSLIB_VERSION)'
	@printf 'LMDB_VERSION=%s\n' '$(LMDB_VERSION)'
	@printf 'HTSLIB_PREFIX=%s\n' '$(HTSLIB_PREFIX)'
	@printf 'LMDB_PREFIX=%s\n' '$(LMDB_PREFIX)'

$(DOWNLOAD_DIR):
	mkdir -p $@

$(OBJ_DIR)/src:
	mkdir -p $@

$(HTSLIB_ARCHIVE): | $(DOWNLOAD_DIR)
	@if command -v curl >/dev/null 2>&1; then \
		curl -LfsS '$(HTSLIB_URL)' -o '$@'; \
	elif command -v wget >/dev/null 2>&1; then \
		wget -q -O '$@' '$(HTSLIB_URL)'; \
	else \
		echo 'curl or wget is required to fetch htslib' >&2; \
		exit 1; \
	fi

$(LMDB_ARCHIVE): | $(DOWNLOAD_DIR)
	@if command -v curl >/dev/null 2>&1; then \
		curl -LfsS '$(LMDB_URL)' -o '$@'; \
	elif command -v wget >/dev/null 2>&1; then \
		wget -q -O '$@' '$(LMDB_URL)'; \
	else \
		echo 'curl or wget is required to fetch LMDB' >&2; \
		exit 1; \
	fi

$(HTSLIB_SRC)/README: $(HTSLIB_ARCHIVE)
	mkdir -p $(VENDOR_SRC_DIR)
	tar -xjf '$<' -C $(VENDOR_SRC_DIR)

$(LMDB_SRC)/Makefile: $(LMDB_ARCHIVE)
	mkdir -p $(VENDOR_SRC_DIR)
	tar -xzf '$<' -C $(VENDOR_SRC_DIR)

$(HTSLIB_LIB): $(HTSLIB_SRC)/README
	cd $(HTSLIB_SRC) && ./configure --prefix='$(HTSLIB_PREFIX)' $(HTSLIB_CONFIGURE_FLAGS)
	$(MAKE) -C $(HTSLIB_SRC)
	$(MAKE) -C $(HTSLIB_SRC) install prefix='$(HTSLIB_PREFIX)'

$(LMDB_LIB): $(LMDB_SRC)/Makefile
	$(MAKE) -C $(LMDB_SRC)
	$(MAKE) -C $(LMDB_SRC) prefix='$(LMDB_PREFIX)' XCFLAGS='-DMDB_MAXKEYSIZE=0' install

$(OBJ_DIR)/src/%.o: src/%.c | $(OBJ_DIR)/src $(HTSLIB_LIB) $(LMDB_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -fvisibility=hidden -fPIC -c $< -o $@

$(BUILD_DIR)/libisam.a: $(ISAM_OBJECTS)
	mkdir -p $(BUILD_DIR)
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(MAIN_OBJECT): src/main.c | $(OBJ_DIR)/src $(HTSLIB_LIB) $(LMDB_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -fvisibility=hidden -fPIC -c $< -o $@

$(BUILD_DIR)/isamtools: $(MAIN_OBJECT) $(BUILD_DIR)/libisam.a $(LMDB_LIB) $(HTSLIB_LIB)
	$(CC) $(LDFLAGS) $(EXE_LDFLAGS) -o $@ $(MAIN_OBJECT) $(BUILD_DIR)/libisam.a $(LMDB_LIB) $(HTSLIB_LIB) $(SYS_LIBS)

.PHONY: all deps test clean distclean print-config