# List of supported compilers
include compilers.mk

# Tools for the shell:
# 1- SEPARATOR, 2- PATH_SEPARATOR, 3- TOUCH, 4- COPY, 5- RENAME, 6- DEL, 7- LOOK_DIR, 8- MKDIR, 9- CPDIR, 10- RMDIR
include $(OS).mk

# Add Toolchain to PATH
export PATH := $(realpath $($(COMPILER)_PATH))$(PATH_SEPARATOR)$(PATH)

# GCC compiler
CC := $($(COMPILER)_CC)

# --------------------------------------------------------------------
WARNINGS_BASICS = 	-Wall -Wextra \
					-Wno-unused-result -Winit-self \
					-Wmissing-declarations -Wredundant-decls -Waggregate-return \
					-Wshadow \
					-Wduplicated-cond -Wlogical-op -Wswitch-default \
					-Wcast-qual -Wcast-align \
					-Wunused-macros -Wpointer-arith -Wmissing-include-dirs \
					-Wmultichar -Woverlength-strings

WARNINGS_FORMAT = -Wformat-nonliteral -Wformat-security -Wformat-signedness -Wformat-y2k
WARNINGS_CONVERSION = -Wsign-conversion -Wfloat-conversion -Wfloat-equal -Wdouble-promotion -Wwrite-strings 
WARNINGS_SUGGEST = -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn
WARNINGS_OPTIMIZE = -Winline -Wunsafe-loop-optimizations -Wdisabled-optimization -Wvector-operation-performance
# --------------------------------------------------------------------
SANITIZES_UNDEFINED = -fsanitize=undefined
# --------------------------------------------------------------------
INSTRUMENTS_linux = -fsplit-stack
INSTRUMENTS_windows = 
INSTRUMENTS_INTEL = -mstack-protector-guard=tls -mmitigate-rop
INSTRUMENTS_ARM =

INSTRUMENTS_STACK = -fstack-check -fstack-protector-strong
INSTRUMENTS_SECTION = -fdata-sections -ffunction-sections -Wl,--gc-sections
INSTRUMENTS_SYSTEM = $(INSTRUMENTS_$(OS))
INSTRUMENTS_PROCESSOR = $(INSTRUMENTS_$(PROCESSOR))
# --------------------------------------------------------------------
OPTIMIZES_BASICS = -finline -fisolate-erroneous-paths-attribute -fgcse-las -fgcse-sm -ftree-partial-pre -fipa-pta
OPTIMIZES_SCHED = 	-fmodulo-sched -fmodulo-sched-allow-regmoves  -fsched-spec-load  -fsel-sched-pipelining -fsel-sched-pipelining-outer-loops \
					-fsched-pressure -fsched-spec-load-dangerous -fselective-scheduling2 -fsched2-use-superblocks
					
OPTIMIZES_LOOP = 	-funsafe-loop-optimizations -ftree-loop-linear -ftree-loop-im -ftree-loop-ivcanon -ftree-loop-if-convert -ftree-loop-distribution -ftree-loop-distribute-patterns \
					-floop-interchange -floop-strip-mine -floop-block -floop-unroll-and-jam -floop-parallelize-all \
					-funswitch-loops -fira-loop-pressure \
					-fivopts -fpredictive-commoning -fgraphite-identity
		
OPTIMIZES_ADVANCED = -fuse-linker-plugin -ffast-math
### -flto

COVERAGES_GENERATE = -fprofile-dir=$(PDIR) -fprofile-generate -fprofile-correction
COVERAGES_USE = -fprofile-dir=$(PDIR) -fprofile-use -fprofile-correction -freorder-blocks-and-partition -freorder-functions -fvariable-expansion-in-unroller
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# Redefine folders with respect to ROOT
SDIR  := $(patsubst $(dir $(SDIR))%,%,$(ROOT)$(SEPARATOR)$(SDIR))
ODIR  := $(patsubst $(dir $(ODIR))%,%,$(ROOT)$(SEPARATOR)$(ODIR))
BDIR  := $(patsubst $(dir $(BDIR))%,%,$(ROOT)$(SEPARATOR)$(BDIR))
BIDIR := $(BDIR)$(SEPARATOR)$(IDIR)
LDIR  := $(patsubst $(dir $(LDIR))%,%,$(ROOT)$(SEPARATOR)$(LDIR))
IDIR  := $(patsubst $(dir $(IDIR))%,%,$(ROOT)$(SEPARATOR)$(IDIR))
PDIR  := $(patsubst $(dir $(PDIR))%,%,$(ROOT)$(SEPARATOR)$(PDIR))
RDIR  := $(patsubst $(dir $(RDIR))%,%,$(ROOT)$(SEPARATOR)$(RDIR))


# Packages
PACKAGES := $(patsubst $(abspath $(SDIR))/%,%,$(realpath $(call LOOK_DIR,$(SDIR))))
PACKAGES := $(filter-out %/$(SDIR),$(PACKAGES))

# source files (including those in packages)
SRCS := $(wildcard $(SDIR)/*.cpp) $(foreach package,$(PACKAGES),$(wildcard $(SDIR)/$(package)/*.cpp))

# object files (including those in packages) 
OBJS := $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCS))
OBJS_LIB := $(filter-out %$(LINK_TARGET).o,$(OBJS))

# Executable file
TARGET := $(BDIR)/$(LINK_TARGET)


# Library Name
LIBRARY_NAME := $(BDIR)/lib$(LIBRARY_NAME)
STATIC_LIBRARY := $(LIBRARY_NAME).a
SHARED_LIBRARY := $(LIBRARY_NAME).so


# Tree to create in init
TREE := $(BDIR) $(ODIR) $(patsubst %,$(ODIR)$(SEPARATOR)%,$(PACKAGES)) $(PDIR)

# Tree for test
TREE_TEST := $(TDIR) $(TDIR)$(SEPARATOR)$(SDIR) $(TDIR)$(SEPARATOR)$(LDIR) $(TDIR)$(SEPARATOR)$(RDIR) $(TDIR)$(SEPARATOR)$(IDIR) $(TDIR)$(SEPARATOR)$(DDIR)

# Tree for temporary tests
TREE_TMPS := 
# ----------------------------------------------------------------------

# --------------------------------------------------------------------
WARNINGS_FLAGS := $(foreach flags,$(WARNINGS),$(WARNINGS_$(flags)))
DEBUGS_FLAGS := -g
SANITIZES_FLAGS := $(WARNINGS_FLAGS) $(DEBUGS_FLAGS) $(foreach flags,$(SANITIZES),$(SANITIZES_$(flags)))
INSTRUMENTS_FLAGS := $(WARNINGS_FLAGS) $(foreach flags,$(INSTRUMENTS),$(INSTRUMENTS_$(flags)))
OPTIMIZES_FLAGS := $(INSTRUMENTS_FLAGS) $(OPTIMIZERS) $(foreach flags,$(OPTIMIZES),$(OPTIMIZES_$(flags)))
RELEASES_FLAGS := $(OPTIMIZES_FLAGS) $(RELEASES)


# Compiler options
CFLAGS_INC := -iquote $(SDIR) -iquote $(IDIR) $(patsubst %,-I%,$(HEADERS_PATH))
CFLAGS = $(CFLAGS_COMMON) $(CFLAGS_INC) $($(MODE)S_FLAGS)

# Linker options
LIBRARIES := $(LIBRARIES_STATIC) $(LIBRARIES_SHARED)
LFLAGS = $(LFLAGS_COMMON) -L$(LDIR) $(patsubst %,-l%,$(LIBRARIES)) $(patsubst lib%,-l%,$(basename $(notdir $(wildcard $(LDIR)/lib*))))
## TO EDIT -Wl,-rpath=$(LDIR)


# ---------------------------- DEPENDENCIES ----------------------------
# Temporary file
TMP_FILE := make_include.mk

# 1- create or overwrite TMP_FILE
TMP_FILE := $(UDIR)$(SEPARATOR)$(TMP_FILE)
TREE_TMP += $(TMP_FILE)
$(shell $(TOUCH) $(TMP_FILE))

# 2- create the dependencies using gcc -MM 
$(foreach file,$(SRCS),$(shell $(CC) -MM -MT $(patsubst $(SDIR)%.cpp,$(ODIR)%.o,$(file)) $(CFLAGS_INC) $(file) >> $(TMP_FILE)))

# include the generated rules
include $(TMP_FILE)
# ----------------------------------------------------------------------




# ---------------------- COMPILING: DO NOT EDIT -----------------------
# static library
libstatic: $(OBJS)
	$(eval $(call MKDIR,$(BIDIR)))
	@gcc-ar rcs $(STATIC_LIBRARY) $(OBJS_LIB)
	$(foreach header,$(LIBRARY_INTERFACE),$(shell $(COPY) $(SDIR)$(SEPARATOR)$(header) $(BIDIR)))
	@echo build static library done

# shared library
libshare: $(OBJS)
	$(eval $(call MKDIR,$(BIDIR)))
	@$(CC) $(CFLAGS) $(COVERAGES_$(COVERAGES)) -o $(SHARED_LIBRARY) $(OBJS_LIB) -shared $(LFLAGS)
	$(foreach header,$(LIBRARY_INTERFACE),$(shell $(COPY) $(SDIR)$(SEPARATOR)$(header) $(BIDIR)))
	@echo build shared library done

# executable
exe: $(OBJS)
	@$(CC) $(CFLAGS) $(COVERAGES_$(COVERAGES)) -o $(TARGET) $(OBJS) $(LFLAGS)
	@echo build executable done

# syntax
syntax:
	@$(CC) -fsyntax-only $(CFLAGS_INC) $(SRCS) 
	@echo syntax OK	

# Compile .cpp files
$(OBJS):
	@echo compiling $@
	@$(CC) -c -o $@ $(patsubst $(ODIR)%.o,$(SDIR)%.cpp,$@) $(CFLAGS) $(COVERAGES_$(COVERAGES))



# Create folders
init: clean_tmps
	$(foreach folder,$(TREE),$(eval $(call MKDIR,$(folder))))
	@echo done init

# Delete TMP files
clean_tmps:
	$(foreach file,$(TREE_TMP),$(shell $(DEL) $(file)))
	
# Delete Coverage profils
clean_coverage:
	$(shell $(RMDIR) $(PDIR))
	@echo done profile cleaning

# Delete objs
clean_objs:
	$(shell $(RMDIR) $(ODIR))

# Delete folders
clean: clean_objs
	$(shell $(RMDIR) $(BDIR))
	@echo done cleaning
# -----------------------------------------------------------------------