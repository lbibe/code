# TO EDIT
NFLLIB_HOME :=
LIBS_HOME   :=



# Project tree
ROOT := .
SDIR := src
UDIR := build
ODIR := obj
PDIR := profiling
IDIR := include
LDIR := lib
BDIR := bin
RDIR := res
DDIR := data

# --------------------------------------------------------------------
# Processors: INTEL or ARM
PROCESSOR = INTEL

# Operating systems: linux or win32
OS = linux

# Compiler flavors: see build/compilers.mk: NATIVE, BWIN
COMPILER := BWIN

# Compiler modes: WARNING or DEBUG or SANITIZE or INSTRUMENT or OPTIMIZE or RELEASE
MODE := RELEASE

# Compiler outputs: syntax or exe or libshare or libstatic
OUTPUT := exe
# --------------------------------------------------------------------
# Target Name : Executable name (and the file containing the main function)
LINK_TARGET := main

# only for Test
TEST_TARGET := $(LINK_TARGET)_test

# Library Name
LIBRARY_NAME := lattice_ibe
LIBRARY_INTERFACE := lattice_ibe.h

# --------------------------------------------------------------------
# Fine-tune Warnings: BASICS CONVERSION SUGGEST OPTIMIZE
WARNINGS =

# Fine-tune Sanitizations: UNDEFINED (to be completed)
SANITIZES = UNDEFINED

# Fine-tune Instruments: STACK SECTION SYSTEM PROCESSOR
INSTRUMENTS = SECTION

# Fine-tune Optimizations: BASICS SCHED LOOP ADVANCED
OPTIMIZES = BASICS SCHED LOOP ADVANCED

# Fine-tune Coverages: GENERATE USE
COVERAGES =
##GENERATE USE

# Fine-tune Release:
RELEASES = -s
# --------------------------------------------------------------------




# --------------------------------------------------------------------
# Headers to be included
HEADERS_PATH =

# Libraries to be included (not mention the ones, *.so and *.a, in the lib folder)
LIBRARIES_SHARED :=
LIBRARIES_STATIC :=

# General compiler options to be fine-tuned
CFLAGS_COMMON := -std=c++14 -I$(NFLLIB_HOME)/include
#-fno-plt -fno-jump-tables -fpic -Wno-missing-field-initializers
## -fpie -fstack-reuse=none

# General linker options to be fine-tuned
LFLAGS_COMMON := -L$(NFLLIB_HOME)/lib -L$(LIBS_HOME) -Wl,-rpath=$(LIBS_HOME):$(NFLLIB_HOME)/lib -lmpfr -lnfllib -lgmp -lpthread
##-static-libubsan

# Optimization options related to the processor
OPTIMIZERS := -O3 -march=native
#-O3 -march=native
##-fsched-stalled-insns=10 -fsched-stalled-insns-dep=49
##-march=skylake -msse2 -mfpmath=sse -malign-data=cacheline -mvzeroupper -maccumulate-outgoing-args

# -------------------------------------------------------------------
default: clean init $(OUTPUT)

build: clean init optimize_cover

# Launch the target
run: clean_tmps
	@./$(TARGET)

# -------------------------------------------------------------------
# Optimize using covering
optimize_cover:
	@$(MAKE) --silent COVERAGES=GENERATE exe run
	@$(MAKE) --silent COVERAGES=USE clean init $(OUTPUT) clean_coverage
# -------------------------------------------------------------------

# include default makefile
include makefile.mk
## Targets: clean clean_objs clean_coverage clean_tmps init init_test exe libstatic libshare libshare_test