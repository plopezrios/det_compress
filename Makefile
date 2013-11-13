#--------------------------------------------------#
# Top-level Makefile for the DET_COMPRESS utility. #
#                                                  #
# PLR 11.2013                                      #
#--------------------------------------------------#

# Basic stuff
SRCDIR  = $(PWD)/src
.SILENT:

# Compilers and options.
FC      = gfortran
FFLAGS  = -O3
CC      = gcc
CFLAGS  = -O3

# Main target
.PHONY: default
default:
	@cd $(SRCDIR) && $(MAKE)

# 'Clean' target
.PHONY: clean
clean:
	@cd $(SRCDIR) && $(MAKE) clean
