#----------------------------------------#
# Makefile for the det_compress utility. #
#                                        #
# PLR 11.2013                            #
#----------------------------------------#

# Basic stuff.
UTILEXE = det_compress
SRCDIR  = $(PWD)
LIBDIR  = $(SRCDIR)/../lib
BINDIR  = $(SRCDIR)/../bin
OBJDIR  = $(SRCDIR)/obj
.SILENT:

# How to include / link the LPSOLVE library.
LPSOLVE_INCDIR  = $(LIBDIR)/lpsolve
LPSOLVE_LIBNAME = lpsolve55
LPSOLVE_LIBDIR  = $(LIBDIR)
LPSOLVE         = $(LPSOLVE_LIBDIR)/lib$(LPSOLVE_LIBNAME).a
LDLPSOLVE_DEPS  = -ldl
INCLUDE_LPSOLVE = -I$(LPSOLVE_INCDIR)

# Compilers and options.
FC      = gfortran
FFLAGS  = -O3
CC      = gcc
CFLAGS  = -O3

# Linker options.
LDLIBS  = -L$(LPSOLVE_LIBDIR) -l$(LPSOLVE_LIBNAME) $(LDLPSOLVE_DEPS)

# Main target.
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# 'Clean' target.
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(BINDIR) ;\
        cd $(LIBDIR)/lpsolve && $(MAKE) clean

# Object-file list
OBJ = \
$(OBJDIR)/casl.o \
$(OBJDIR)/det_compress.o \
$(OBJDIR)/tools.o
# End object-file list

# Module dependencies
$(OBJDIR)/casl.o : $(SRCDIR)/casl.f90
casl.o_name = casl

DET_COMPRESS_DEP = \
$(OBJDIR)/casl.o \
$(OBJDIR)/tools.o
$(OBJDIR)/det_compress.o : $(SRCDIR)/det_compress.f90 $(DET_COMPRESS_DEP)
det_compress.o_name = det_compress

$(OBJDIR)/tools.o : $(SRCDIR)/tools.f90
tools.o_name = tools
# End module dependencies

# How to make the objects.
$(LPSOLVE):
	@cd $(LIBDIR)/lpsolve && $(MAKE)
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        echo "FC $($(@F)_name).f90" ;\
        cd $(OBJDIR) ;\
        $(FC) -c $(FFLAGS) $(SRCDIR)/$($(@F)_name).f90
$(OBJDIR)/lp2fortran.o: $(SRCDIR)/lp2fortran.c
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        echo "CC lp2fortran.c" ;\
        cd $(OBJDIR) ;\
        $(CC) $(CFLAGS) $(CFLAGS_FC_INTERFACE) -c $(SRCDIR)/lp2fortran.c\
        $(INCLUDE_LPSOLVE)
$(BINDIR)/$(UTILEXE): $(LPSOLVE) $(OBJDIR)/lp2fortran.o $(OBJ)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        echo "FC $(UTILEXE)" ;\
        $(FC) $(FFLAGS) -o $(@) $(OBJ) $(OBJDIR)/lp2fortran.o $(LDLIBS)
