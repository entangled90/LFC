################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq rmpmk
.PHONY: all


# main programs and required modules 

MAIN = harmonicMetropolis 
# importanza harmonic_gamma buffon

RANDOM = ranlxs ranlxd gauss pdf

START = start 

EXTRAS = varie my_struct

DERIVATIVE =
INTEGRATION = integration montecarlo
HARMONIC = harmonic

MODULES = $(RANDOM) $(START) $(EXTRAS) $(DERIVATIVE) $(INTEGRATION) $(HARMONIC)


# search path for modules

MDIR = ../modules

VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/extras:$(MDIR)/derivative:$(MDIR)/integration:$(MDIR)/harmonic



# additional include directories

INCPATH = ../include


# additional libraries to be included 
 
LIBS = m

LIBPATH =


# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -Wall -O -g  -std=gnu99 -pedantic -DSSE2
 

############################## do not change ###################################

SHELL=/bin/bash

CC=gcc

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(CC) -MM -ansi $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"

rmpmk: 
	@ -rm -rf *.d *.o *.tmp

# clean directory 

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
