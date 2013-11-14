#--------------------------------------------------#
# Top-level Makefile for the DET_COMPRESS utility. #
#                                                  #
# PLR 11.2013                                      #
#--------------------------------------------------#

# Basic stuff
SRCDIR   = $(PWD)/src
WDEMODIR = $(PWD)/interface/write_mdet/src
RDEMODIR = $(PWD)/interface/read_cmdet/src
.SILENT:

# Compilers and options.
FC      = gfortran
FFLAGS  = -O3
CC      = gcc
CFLAGS  = -O3

# Main target
.PHONY: default
default: src read_demo write_demo

# Compilation targets.
.PHONY: src
src:
	@cd $(SRCDIR) && $(MAKE)

.PHONY: read_demo
read_demo:
	@cd $(RDEMODIR) && $(MAKE)

.PHONY: write_demo
write_demo:
	@cd $(WDEMODIR) && $(MAKE)

# 'Clean' targets.
.PHONY: clean
clean: clean_src clean_read_demo clean_write_demo

.PHONY: clean_src
clean_src:
	@cd $(SRCDIR) && $(MAKE) clean

.PHONY: clean_read_demo
clean_read_demo:
	@cd $(RDEMODIR) && $(MAKE) clean

.PHONY: clean_write_demo
clean_write_demo:
	@cd $(WDEMODIR) && $(MAKE) clean
