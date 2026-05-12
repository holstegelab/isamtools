# Migrating an htslib package to libhtsisam

This guide walks through making an existing htslib-based program transparently
read `.isam` files alongside regular SAM/BAM/CRAM, with no behavioural change
on plain BAM/CRAM input.

The compat layer (`include/htsisam/sam.h` + `htsisam_dispatch.c`, compiled into
`libisam.a`) intercepts a small set of htslib entry points (`hts_open`,
`sam_read1`, `sam_index_load`, `sam_itr_querys/queryi/next`, etc.). At runtime
it picks between real htslib and libisam based on the path suffix, using a
tagged-pointer scheme so the cost on the BAM/CRAM path is one bit-test per call.

A worked example for the `otter` package lives at
[`api_compat_example/Makefile.htsisam`](../api_compat_example/Makefile.htsisam).

## What gets intercepted

Only these symbols are redirected:

```
hts_open  hts_open_format  sam_open  sam_close  hts_close
sam_hdr_read  sam_hdr_write
sam_read1  sam_write1
sam_index_load  sam_index_load2  sam_index_load3  hts_idx_destroy
sam_itr_querys  sam_itr_queryi  sam_itr_next  sam_itr_multi_next
hts_itr_destroy
bam_itr_querys  bam_itr_queryi  bam_itr_next
```

Everything else from `<htslib/sam.h>` (record helpers like `bam_get_seq`,
`bam_aux_get`, `bam_endpos`, `bam_init1`, header helpers like
`sam_hdr_line_name`, faidx, kstring, BGZF) is unchanged. Records are
real `bam1_t`, headers are real `sam_hdr_t`.

## Choosing a migration path

Pick the lowest level that still works for your project.

### Path A — compile-time `-include`, zero source edits

Best when the package builds with a Makefile/CMake you control and you don't
want to touch sources at all.

Add one compile flag and three link inputs.

**Compile** (every translation unit):

```
-I/path/to/isamtools/include   \
-include htsisam/sam.h
```

The `-include` flag pre-includes our header into every TU before any
user-supplied `<htslib/sam.h>`. Inside our header, after `<htslib/sam.h>` is
included, `sam_open`, `sam_read1`, etc. are `#define`d to the wrapping
`htsisam_*` symbols.

**Link**:

```
/path/to/isamtools/build/libisam.a       \
/path/to/isamtools/build/vendor/htslib/lib/libhts.a   \
/path/to/isamtools/build/vendor/lmdb/lib/liblmdb.a    \
/path/to/isamtools/build/vendor/bzip2/lib/libbz2.a    \
-lz -llzma -lpthread -lm
```

The vendored static archives in the isamtools build are not built with
`-fPIC`, so add `-no-pie` to the executable link line. (libisam.a itself is
PIC, so it can also be embedded into `-fPIC` shared libraries if you prefer to
build a `libhtsisam.so` and dynamically link against that.)

That is the entire migration. No source files changed.

### Path B — change one include in your project

Best when one header in your project includes htslib and pulls it into many
TUs. Replace there:

```c
/* before */
#include <htslib/sam.h>

/* after */
#include <htsisam/sam.h>
```

The htsisam header itself includes `<htslib/sam.h>`, so symbols like `bam1_t`
remain available. Build and link as in Path A but you can drop the `-include`
flag.

### Path C — explicit `htsisam_*` calls, no macro shadow

Best when:

- You want `htsisam_*` to show up in stack traces and tooling.
- Your codebase has macros named `sam_open` (or similar) that would collide
  with the dropin macros.
- You want to call htslib and libisam side-by-side in the same TU on the same
  pointer types.

Define `HTSISAM_NO_DROPIN` before including the header and call the wrappers
directly:

```c
#define HTSISAM_NO_DROPIN
#include <htsisam/sam.h>
...
htsFile *fp = htsisam_hts_open(path, "r");
sam_hdr_t *hdr = htsisam_sam_hdr_read(fp);
while (htsisam_sam_read1(fp, hdr, b) >= 0) { ... }
htsisam_sam_close(fp);
```

Same link inputs as Path A.

## Worked example: otter (Path A)

Otter's upstream Makefile compiles `src/*.cpp` and links `-lhts`. The migrated
build adds three things:

```make
ISAM_DIR     ?= /path/to/isamtools
ISAM_BUILD   ?= $(ISAM_DIR)/build

HTSISAM_CPPFLAGS = -I$(ISAM_DIR)/include \
                   -I$(ISAM_BUILD)/vendor/htslib/include \
                   -include htsisam/sam.h

EXE_LDFLAGS = -no-pie

LINK_LIBS   = $(ISAM_BUILD)/libisam.a                          \
              $(ISAM_BUILD)/vendor/htslib/lib/libhts.a         \
              $(ISAM_BUILD)/vendor/lmdb/lib/liblmdb.a          \
              $(ISAM_BUILD)/vendor/bzip2/lib/libbz2.a          \
              -lwfacpp -lz -llzma -lpthread -lm -fopenmp

$(ODIR)/%.o: $(SDIR)/%.cpp | $(ODIR)
	$(CXX) $(CXX_FLAGS) $(HTSISAM_CPPFLAGS) -MMD -MP -c -o $@ $< $(INC)

$(ODIR)/$(OUT): $(OBJS) | $(ODIR)
	$(CXX) $(CXX_FLAGS) $(EXE_LDFLAGS) -o $@ $^ $(INC) $(LINK_LIBS)
```

Verification: run otter on the same data through both inputs and check
md5-equality of the output.

```sh
otter assemble -b regions.bed sample.bam   | md5sum
otter assemble -b regions.bed sample.isam  | md5sum
```

Both produce identical bytes.

## Caveats and limits

**Write mode is read-only for `.isam`.** `sam_open("path.isam", "wb")` is
refused with a clear error. Writing a `.isam` is an import operation
(`isamtools import`), not a stream of `sam_write1` calls. Plain BAM/CRAM write
modes pass through unchanged.

**Don't reach into struct internals on `.isam` paths.** Code that touches
`fp->fp.bgzf`, `fp->is_bgzf`, `fp->fn`, or other `htsFile` fields will misread
on `.isam` (the pointer is tagged libisam state, not a real `htsFile`). The
common public API (`sam_*`, `bam_*` helpers) is fine. If your code does this,
gate it on a path check:

```c
if (!htsisam_is_isam_path(path)) {
    /* htslib-internal access here */
}
```

**Path detection is by suffix** (`*.isam`, case-insensitive). A magic-byte
sniff fallback is a small follow-up if you have ambiguous filenames.

**Iterator pointers are tagged.** `hts_itr_t *` returned from `sam_itr_querys`
on a `.isam` is not a real htslib iterator. Don't pass it to
`hts_itr_next` directly, only via `sam_itr_next` (which is what the macro
shims in `<htsisam/sam.h>` redirect to).

**Threads.** `hts_set_threads(fp, n)` on a `.isam` is a no-op today; threading
on the canonical-store side is exposed via `isamOptions.n_threads` to programs
that use the native libisam API. Worth wiring through the compat layer if your
consumer relies on it.

**`sam_itr_next` return convention.** htslib returns the record's serialised
size on success (positive int); some consumers test `> 0`, others `>= 0`. The
compat layer matches htslib (positive on success) so both idioms work.

## Migration checklist

Before declaring done:

- [ ] Build the project with `-include htsisam/sam.h` (or with the include
      swap from Path B). Compilation succeeds with no `redefined` errors;
      `sam_open` / `sam_close` warnings about prior macros are expected and
      benign.
- [ ] Verify `nm <binary> | grep htsisam_` shows the wrappers linked in.
- [ ] Run the unmodified test suite against a plain BAM. Output must be
      bit-identical to the pre-migration build.
- [ ] Run the same test against the `.isam` of the same data. Output must
      match the BAM run modulo any deliberate filtering differences (e.g.
      records that aren't in the canonical store and were `--allow-missing`
      during import).
- [ ] If your tool writes BAM/CRAM, confirm the write path still works
      (it goes through real htslib unchanged).

## Building libisam.a for distribution

From the isamtools project root:

```sh
make                              # builds build/libisam.a + vendored deps
ls build/libisam.a                # the static archive consumers link
ls build/vendor/htslib/lib/libhts.a build/vendor/lmdb/lib/liblmdb.a
ls include/htsisam/sam.h          # the dropin header consumers `-include`
```

Distribute `libisam.a`, the htslib + lmdb static archives, and the
`include/htsisam` and `include/isam` directories. Consumers add the include
directory to `-I` and the three archives to their link line.

## When the wrapper is not enough

Two cases push beyond the dropin layer:

1. **Programs that interpose into htslib internals** — samtools' own internal
   subcommands, packages that read `htsFile->fn` or call `hts_get_format()`.
   Those need code-level attention or a shim that synthesises a fake
   `htsFile` with sensible field values, which the current compat layer does
   not do.

2. **Programs that need `.isam`-aware features** — annotation merging,
   subset-filtered iteration, multi-view joining. The htslib API doesn't
   express those concepts; use the native libisam C API
   (`include/isam/isam.h`, `isam_open`, `isam_itr_query_subset`, etc.)
   directly for those code paths.

For most read-path consumers (variant callers, coverage tools, assemblers,
visualisation backends), Path A is enough.
