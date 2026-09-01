#===================================================================
#
#             Makefile for Vampire Data Converter (vdc)
#
#===================================================================

# Compilers
GCC=g++

# LIBS
LIBS=-lstdc++

# Flags
#GCC_CFLAGS=-O3 -std=c++0x
GCC_CFLAGS=-O3 -std=c++0x -fopenmp -funroll-loops #-mavx2

#-------------------------------------------------------------------
# Optional FFTW3 for --stray-field (2D magnetostatic FFT).
#
# Auto-detects a working fftw3 header+library. If FFTW is missing,
# vdc still compiles; --stray-field then prints an error telling you
# to recompile with FFTW.
#
#   make              # auto
#   make FFTW=0       # force off
#   make FFTW=1       # force on (fails at compile/link if missing)
#-------------------------------------------------------------------
FFTW ?= auto

FFTW_PKG_CFLAGS := $(shell pkg-config --cflags fftw3 2>/dev/null)
FFTW_PKG_LIBS   := $(shell pkg-config --libs fftw3 2>/dev/null)
ifeq ($(FFTW_PKG_LIBS),)
  FFTW_PKG_LIBS := -lfftw3
endif

ifeq ($(filter clean purge,$(MAKECMDGOALS)),)
ifeq ($(FFTW),0)
  FFTW_CFLAGS :=
  FFTW_LIBS :=
  FFTW_STATUS := disabled
else ifeq ($(FFTW),1)
  FFTW_CFLAGS := -DVDC_FFTW $(FFTW_PKG_CFLAGS)
  FFTW_LIBS := $(FFTW_PKG_LIBS)
  FFTW_STATUS := forced-on
else
  FFTW_PROBE := $(shell mkdir -p obj && echo 'int main(){ fftw_free(fftw_malloc(8)); return 0; }' | $(GCC) $(FFTW_PKG_CFLAGS) -include fftw3.h -x c++ - $(FFTW_PKG_LIBS) -o obj/.fftw_probe >/dev/null 2>&1 && echo yes)
  ifeq ($(FFTW_PROBE),yes)
    FFTW_CFLAGS := -DVDC_FFTW $(FFTW_PKG_CFLAGS)
    FFTW_LIBS := $(FFTW_PKG_LIBS)
    FFTW_STATUS := enabled
  else
    FFTW_CFLAGS :=
    FFTW_LIBS :=
    FFTW_STATUS := not-found
  endif
endif
else
  FFTW_CFLAGS :=
  FFTW_LIBS :=
  FFTW_STATUS := skipped
endif

# Objects
OBJECTS= \
obj/binary.o \
obj/cells.o \
obj/colour.o \
obj/colourmaps.o \
obj/command.o \
obj/coords.o \
obj/data.o \
obj/functions.o \
obj/main.o \
obj/non_magnetic.o \
obj/povray.o \
obj/povray-cells.o \
obj/read.o \
obj/spins.o \
obj/ssc.o \
obj/stray_field.o \
obj/text.o \
obj/vtk.o \
obj/xyz.o

EXECUTABLE=vdc

all: $(OBJECTS) gcc

GCCDB_OBJECTS=$(OBJECTS:.o=_gdb.o)

GCC_DBCFLAGS= -g -Wall -Wextra -O0 -fbounds-check -pedantic -std=c++11 -Wno-long-long
GCC_DBLFLAGS= -g -std=c++11 -fbounds-check

OMP_FLAGS= -fopenmp

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(OBJECTS) $(GCC_CFLAGS) $(LIBS) $(FFTW_LIBS) -o $(EXECUTABLE)
	@echo "vdc: FFTW3 $(FFTW_STATUS)"

gcc-debug: $(GCCDB_OBJECTS)
	$(GCC) $(GCC_DBLFLAGS) $(OMP_FLAGS) $(GCCDB_OBJECTS) $(FFTW_LIBS) -o $(EXECUTABLE)-debug
	@echo "vdc-debug: FFTW3 $(FFTW_STATUS)"

$(GCCDB_OBJECTS): obj/%_gdb.o: ./%.cpp
	$(GCC) -c -o $@ $(GCC_DBCFLAGS) $(OMP_FLAGS) $(FFTW_CFLAGS) $<

$(OBJECTS): obj/%.o: ./%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $(FFTW_CFLAGS) $<

clean:
	@rm -f obj/*.o obj/.fftw_probe

purge:
	@rm -f obj/*.o obj/.fftw_probe
	@rm -f $(EXECUTABLE) $(EXECUTABLE)-debug
