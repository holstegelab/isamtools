.DEFAULT_GOAL := all

CC ?= cc
AR ?= ar
RANLIB ?= ranlib
MAKE ?= make
BUILD_DIR ?= build

HTSLIB_VERSION ?= 1.23.1
LMDB_VERSION ?= LMDB_0.9.31
BZIP2_VERSION ?= 1.0.8

DOWNLOAD_DIR := $(BUILD_DIR)/downloads
VENDOR_SRC_DIR := $(BUILD_DIR)/vendor/src
OBJ_DIR := $(BUILD_DIR)/obj

HTSLIB_URL ?= https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2
LMDB_URL ?= https://github.com/LMDB/lmdb/archive/$(LMDB_VERSION).tar.gz
BZIP2_URL ?= https://sourceware.org/pub/bzip2/bzip2-$(BZIP2_VERSION).tar.gz

HTSLIB_ARCHIVE := $(DOWNLOAD_DIR)/htslib-$(HTSLIB_VERSION).tar.bz2
LMDB_ARCHIVE := $(DOWNLOAD_DIR)/lmdb-$(LMDB_VERSION).tar.gz
BZIP2_ARCHIVE := $(DOWNLOAD_DIR)/bzip2-$(BZIP2_VERSION).tar.gz

HTSLIB_SRC := $(VENDOR_SRC_DIR)/htslib-$(HTSLIB_VERSION)
LMDB_SRC := $(VENDOR_SRC_DIR)/lmdb-$(LMDB_VERSION)/libraries/liblmdb
BZIP2_SRC := $(VENDOR_SRC_DIR)/bzip2-$(BZIP2_VERSION)

HTSLIB_PREFIX := $(abspath $(BUILD_DIR)/vendor/htslib)
LMDB_PREFIX := $(abspath $(BUILD_DIR)/vendor/lmdb)
BZIP2_PREFIX := $(abspath $(BUILD_DIR)/vendor/bzip2)

HTSLIB_LIB := $(HTSLIB_PREFIX)/lib/libhts.a
LMDB_LIB := $(LMDB_PREFIX)/lib/liblmdb.a
BZIP2_LIB := $(BZIP2_PREFIX)/lib/libbz2.a

HTSLIB_CONFIGURE_FLAGS ?= --disable-libcurl --disable-plugins --disable-s3 --disable-gcs --without-libdeflate
HTSLIB_CONFIGURE_CPPFLAGS ?= -I$(BZIP2_PREFIX)/include
HTSLIB_CONFIGURE_LDFLAGS ?= -L$(BZIP2_PREFIX)/lib
SYS_LIBS ?= -lz -llzma -lm -lpthread -lzstd
EXE_LDFLAGS ?= -no-pie

WARNFLAGS := -Wall -Wextra -Wno-unused-parameter -Wshadow
CFLAGS ?= -O3 -g
CPPFLAGS += -Iinclude -Isrc -I$(HTSLIB_PREFIX)/include -I$(LMDB_PREFIX)/include \
	-D_POSIX_C_SOURCE=200809L -D_GNU_SOURCE
LDFLAGS +=

ISAM_SOURCES := \
	src/common.c \
	src/aux_tags.c \
	src/backend_lmdb.c \
	src/hash.c \
	src/rebuild.c \
	src/store.c \
	src/view.c \
	src/cmd_index.c \
	src/cmd_import.c \
	src/cmd_check.c \
	src/cmd_sort.c \
	src/cmd_view.c \
	src/cmd_annotate.c \
	src/cmd_store_stats.c \
	src/cmd_audit_tags.c \
	src/cmd_update.c \
	src/cmd_rehint.c \
	src/api.c \
	src/htsisam_dispatch.c

ISAM_OBJECTS := $(patsubst src/%.c,$(OBJ_DIR)/src/%.o,$(ISAM_SOURCES))
MAIN_OBJECT := $(OBJ_DIR)/src/main.o

all: $(BUILD_DIR)/libisam.a $(BUILD_DIR)/isamtools

deps: $(HTSLIB_LIB) $(LMDB_LIB) $(BZIP2_LIB)

packages: deps

test: all
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/functional.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/lmdb_store.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/lmdb_alignment.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) bash tests/tag_audit.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) bash tests/tag_policy.sh
	TOOL=$(abspath $(BUILD_DIR)/isamtools) tests/lmdb_compat.sh

clean:
	rm -rf $(OBJ_DIR) $(BUILD_DIR)/libisam.a $(BUILD_DIR)/isamtools

distclean:
	rm -rf $(BUILD_DIR)

print-config:
	@printf 'BUILD_DIR=%s\n' '$(BUILD_DIR)'
	@printf 'HTSLIB_VERSION=%s\n' '$(HTSLIB_VERSION)'
	@printf 'LMDB_VERSION=%s\n' '$(LMDB_VERSION)'
	@printf 'BZIP2_VERSION=%s\n' '$(BZIP2_VERSION)'
	@printf 'HTSLIB_PREFIX=%s\n' '$(HTSLIB_PREFIX)'
	@printf 'LMDB_PREFIX=%s\n' '$(LMDB_PREFIX)'
	@printf 'BZIP2_PREFIX=%s\n' '$(BZIP2_PREFIX)'

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

$(BZIP2_ARCHIVE): | $(DOWNLOAD_DIR)
	@if command -v curl >/dev/null 2>&1; then \
		curl -LfsS '$(BZIP2_URL)' -o '$@'; \
	elif command -v wget >/dev/null 2>&1; then \
		wget -q -O '$@' '$(BZIP2_URL)'; \
	else \
		echo 'curl or wget is required to fetch bzip2' >&2; \
		exit 1; \
	fi

$(HTSLIB_SRC)/README: $(HTSLIB_ARCHIVE)
	mkdir -p $(VENDOR_SRC_DIR)
	tar -xjf '$<' -C $(VENDOR_SRC_DIR)
	tar -xOf $< htslib-$(HTSLIB_VERSION)/htscodecs_bundled.mk > $(HTSLIB_SRC)/htscodecs_bundled.mk

$(LMDB_SRC)/Makefile: $(LMDB_ARCHIVE)
	mkdir -p $(VENDOR_SRC_DIR)
	tar -xzf '$<' -C $(VENDOR_SRC_DIR)

$(BZIP2_SRC)/Makefile: $(BZIP2_ARCHIVE)
	mkdir -p $(VENDOR_SRC_DIR)
	tar -xzf '$<' -C $(VENDOR_SRC_DIR)

$(HTSLIB_LIB): $(HTSLIB_SRC)/README $(BZIP2_LIB)
	cd $(HTSLIB_SRC) && CPPFLAGS='$(HTSLIB_CONFIGURE_CPPFLAGS)' LDFLAGS='$(HTSLIB_CONFIGURE_LDFLAGS)' ./configure --prefix='$(HTSLIB_PREFIX)' $(HTSLIB_CONFIGURE_FLAGS)
	rm -f $(HTSLIB_SRC)/htscodecs.mk
	tar -xOf $(HTSLIB_ARCHIVE) htslib-$(HTSLIB_VERSION)/htscodecs_bundled.mk > $(HTSLIB_SRC)/htscodecs_bundled.mk
	cp $(HTSLIB_SRC)/htscodecs_bundled.mk $(HTSLIB_SRC)/htscodecs.mk
	$(MAKE) -C $(HTSLIB_SRC) CPPFLAGS='$(HTSLIB_CONFIGURE_CPPFLAGS)' LDFLAGS='$(HTSLIB_CONFIGURE_LDFLAGS)' lib-static
	mkdir -p $(HTSLIB_PREFIX)/include/htslib $(HTSLIB_PREFIX)/lib
	cp $(HTSLIB_SRC)/htslib/*.h $(HTSLIB_PREFIX)/include/htslib/
	cp $(HTSLIB_SRC)/libhts.a $(HTSLIB_PREFIX)/lib/

$(LMDB_LIB): $(LMDB_SRC)/Makefile
	$(MAKE) -C $(LMDB_SRC)
	$(MAKE) -C $(LMDB_SRC) prefix='$(LMDB_PREFIX)' XCFLAGS='-DMDB_MAXKEYSIZE=0' install

$(BZIP2_LIB): $(BZIP2_SRC)/Makefile
	$(MAKE) -C $(BZIP2_SRC) clean
	$(MAKE) -C $(BZIP2_SRC) CC='$(CC)' AR='$(AR)' RANLIB='$(RANLIB)' CFLAGS='$(CFLAGS) -fPIC'
	mkdir -p $(BZIP2_PREFIX)/include $(BZIP2_PREFIX)/lib
	cp $(BZIP2_SRC)/bzlib.h $(BZIP2_PREFIX)/include/
	cp $(BZIP2_SRC)/libbz2.a $(BZIP2_PREFIX)/lib/

$(OBJ_DIR)/src/%.o: src/%.c | $(OBJ_DIR)/src $(HTSLIB_LIB) $(LMDB_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -fvisibility=hidden -fPIC -c $< -o $@

$(BUILD_DIR)/libisam.a: $(ISAM_OBJECTS)
	mkdir -p $(BUILD_DIR)
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(MAIN_OBJECT): src/main.c | $(OBJ_DIR)/src $(HTSLIB_LIB) $(LMDB_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(WARNFLAGS) -fvisibility=hidden -fPIC -c $< -o $@

$(BUILD_DIR)/isamtools: $(MAIN_OBJECT) $(BUILD_DIR)/libisam.a $(LMDB_LIB) $(HTSLIB_LIB) $(BZIP2_LIB)
	$(CC) $(LDFLAGS) $(EXE_LDFLAGS) -o $@ $(MAIN_OBJECT) $(BUILD_DIR)/libisam.a $(LMDB_LIB) $(HTSLIB_LIB) $(BZIP2_LIB) $(SYS_LIBS)

.PHONY: all deps packages test clean distclean print-config
