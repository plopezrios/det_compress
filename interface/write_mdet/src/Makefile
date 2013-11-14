#-------------------------------------------------#
# Makefile for the write_demo demonstration tool. #
#                                                 #
# PLR 11.2013                                     #
#-------------------------------------------------#

# Basic stuff.
UTILEXE = write_demo
SRCDIR  = $(PWD)
LIBDIR  = $(SRCDIR)/../lib
BINDIR  = $(SRCDIR)/../bin
OBJDIR  = $(SRCDIR)/obj
.SILENT:

# Compilers and options.
FC      = gfortran
FFLAGS  = -O3

# Main target.
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# 'Clean' target.
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(BINDIR)

# Object-file list
OBJ = \
$(OBJDIR)/casl.o \
$(OBJDIR)/write_demo.o
# End object-file list

# Module dependencies
$(OBJDIR)/casl.o : $(SRCDIR)/casl.f90
casl.o_name = casl

READ_DEMO_DEP = \
$(OBJDIR)/casl.o
$(OBJDIR)/write_demo.o : $(SRCDIR)/write_demo.f90 $(DET_COMPRESS_DEP)
write_demo.o_name = write_demo
# End module dependencies

# How to make the objects.
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        echo "FC $($(@F)_name).f90" ;\
        cd $(OBJDIR) ;\
        $(FC) -c $(FFLAGS) $(SRCDIR)/$($(@F)_name).f90
$(BINDIR)/$(UTILEXE): $(OBJ)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        echo "FC $(UTILEXE)" ;\
        $(FC) $(FFLAGS) -o $(@) $(OBJ)
