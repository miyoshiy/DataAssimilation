# Master makefile for Numerical Recipes

CCKR = cc
CCANSI = gcc
#gcc has been found to be the best bet

# Important system parameters
#FFLAGS =
CFLAGSANSI = -I../include
CFLAGSKR = -I../include
# Use -O if you trust your optimizer, and
# -g to make a debugging version of this library
# See the document COMPILING and your system documentation for information
# on values to use here.

# Library names . . .
CKRLIB = librecipes_ckr.a
CANSILIB = librecipes_c.a
FLIB = librecipes_f.a
# and destination.
LIBDIR = /usr/local/lib

TEST = /bin/test
# Alternate location on some systems (or trust the system path) :
# TEST = /usr/ucb/test

AR = /bin/ar
# Alternate location on some systems (or trust the system path) :
# AR = /usr/ccs/bin/ar
ARFLAGS = crv

RANLIB = $(AR) -ts
# IMPORTANT: Ranlib is now obsolete. The above only for compatibility. If
# this creates problems, try following substitution in this make machinery :
# RANLIB = ranlib

MAKE = /bin/make
# Alternate location on some systems (or trust the system path) :
# MAKE = /usr/ccs/bin/make

CBIN = ../bin
FBIN = ../bin
# This is internal

SHELL = /bin/sh
# Just to make sure

COMMON_ENV = "LIBDIR = $(LIBDIR)" \
	"FC = $(FC)" "CCKR = $(CCKR)" "CCANSI = $(CCANSI)" \
	"FFLAGS = $(FFLAGS)" \
	"CFLAGSANSI = $(CFLAGSANSI)" "CFLAGSKR = $(CFLAGSKR)" \
	"FBIN = $(FBIN)" "CBIN = $(CBIN)" "FLIB = $(FLIB)" \
	"CANSILIB = $(CANSILIB)" "CKRLIB = $(CKRLIB)" \
	"RANLIB = $(RANLIB)" \
	"AR = $(AR)" "ARFLAGS = $(ARFLAGS)" "SHELL = $(SHELL)"

MAKE_LIST = make_list

all: distributed libs install

clear_dist:
	/bin/rm -f $(MAKE_LIST)
	cp /dev/null $(MAKE_LIST)

fortran_if_have:
	-if $(TEST) -d ../recipes_f; then \
	echo '../recipes_f' >> $(MAKE_LIST); fi

c_ansi_if_have:
	-if $(TEST) -d ../recipes_c-ansi; then \
	echo '../recipes_c-ansi' >> $(MAKE_LIST); fi

c_kr_if_have:
	-if $(TEST) -d ../recipes_c-kr; then \
	echo '../recipes_c-kr' >> $(MAKE_LIST); fi

distributed: clear_dist fortran_if_have c_ansi_if_have c_kr_if_have

libs:
	echo Making Numerical Recipes
	for distrib in `cat $(MAKE_LIST)`; do \
	for library in $$distrib/recipes; do \
	(cd $$library && $(MAKE) all $(COMMON_ENV)); done; done

install: 
	echo Making Install
	for distrib in `cat $(MAKE_LIST)`; do \
	for library in $$distrib/recipes; do \
	(cd $$library && $(MAKE) install $(COMMON_ENV)); done; done

demobins: 
	echo Making Numerical Recipes Demonstration output
	for distrib in `cat $(MAKE_LIST)`; do \
	for bindir in $$distrib/demo/bin; do \
	(cd $$bindir && $(MAKE) all $(COMMON_ENV)); done; done

results:
	echo Making Numerical Recipes Demonstration output
	for distrib in `cat $(MAKE_LIST)`; do \
	for resdir in $$distrib/demo/results; do \
	(cd $$resdir && $(MAKE) all $(COMMON_ENV)); done; done

clean: 
	echo Making Clean
	for distrib in `cat $(MAKE_LIST)`; do \
	for library in $$distrib/recipes; do \
	(cd $$library && $(MAKE) clean $(COMMON_ENV)); done; done
