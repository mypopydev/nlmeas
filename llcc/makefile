# Copyright 2011 Jose-Luis Lisani <joseluis.lisani@uib.es>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.


# C source code
CSRC	= library/io_png.c library/MCMMondelli/FDS_MCM.c \
          library/gaussianGetreuer/gaussian_conv_sii.c \
          library/gaussianGetreuer/gaussian_conv.c

# C++ source code
CXXSRC	= library/io_RGB.cpp 

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= llcc

default	: $(BIN)

# use DEBUG mode with `make DEBUG=1`
ifdef DEBUG
  # C optimization flags
  COPT = -g
else
  # C optimization flags
  COPT	= -O3 -ftree-vectorize -funroll-loops
endif


# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT)  -std=c99
# C++ compilation flags
# option -DNO_XML needed to compile fast bilateral filter 
# implementation by Sylvain Paris and Frédo Durand
CXXFLAGS	= $(CXXOPT) -ansi -DNO_XML
# link flags
LDFLAGS	= -lpng  


CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas


# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS) 
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $<  $(CXXFLAGS) 

# link all the object code
$(BIN): % : %.o $(OBJ) 
	$(CXX) -o $@ $^ $(LDFLAGS) 

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
distclean	: clean
	$(RM) $(BIN) *.o *~ 

