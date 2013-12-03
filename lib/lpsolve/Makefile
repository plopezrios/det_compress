#--------------------------------------------------#
# Makefile to create a static library for LPSOLVE. #
#                                                  #
# PLR 11.2013                                      #
#--------------------------------------------------#

# Basic stuff
SRCDIR         = $(PWD)
OBJDIR         = $(SRCDIR)/obj
LIBDIR         = $(SRCDIR)/..
LPSOLVE_AR     = liblpsolve55.a
LPSOLVE        = $(LIBDIR)/$(LPSOLVE_AR)
COLAMD_SRCDIR  = $(LIBDIR)/colamd
COLAMD         = $(LIBDIR)/libcolamd.a
INCLUDE_COLAMD = -I$(COLAMD_SRCDIR)
.SILENT:

# Defaults for system-dependent variables
AR     = ar crs
CC     = gcc
CFLAGS = -O3

# Main target
.PHONY: default
default: $(LPSOLVE)

# Target for deleting object files
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(LPSOLVE) ;\
        cd $(COLAMD_SRCDIR) && $(MAKE) clean

# Dependencies.
$(COLAMD):
	@cd $(COLAMD_SRCDIR) && $(MAKE)

# Object-file list
OBJ = \
$(OBJDIR)/lp_MDO.o \
$(OBJDIR)/commonlib.o \
$(OBJDIR)/fortify.o \
$(OBJDIR)/ini.o \
$(OBJDIR)/lp_crash.o \
$(OBJDIR)/lp_Hash.o \
$(OBJDIR)/lp_lib.o \
$(OBJDIR)/lp_LUSOL.o \
$(OBJDIR)/lp_matrix.o \
$(OBJDIR)/lp_mipbb.o \
$(OBJDIR)/lp_MPS.o \
$(OBJDIR)/lp_params.o \
$(OBJDIR)/lp_presolve.o \
$(OBJDIR)/lp_price.o \
$(OBJDIR)/lp_pricePSE.o \
$(OBJDIR)/lp_report.o \
$(OBJDIR)/lp_rlp.o \
$(OBJDIR)/lp_scale.o \
$(OBJDIR)/lp_simplex.o \
$(OBJDIR)/lp_SOS.o \
$(OBJDIR)/lp_utils.o \
$(OBJDIR)/lp_wlp.o \
$(OBJDIR)/lusol.o \
$(OBJDIR)/mmio.o \
$(OBJDIR)/myblas.o \
$(OBJDIR)/yacc_read.o
# End object-file list

# Handle native compilation
BUNCH_OF_DEFINES = -DYY_NEVER_INTERACTIVE -DPARSER_LP\
 -DINVERSE_ACTIVE=INVERSE_LUSOL -DRoleIsExternalInvEngine

# How to compile things
$(LPSOLVE): $(COLAMD) $(OBJ)
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        echo "AR LPSOLVE" ;\
        rm -f $(LPSOLVE) ;\
        $(AR) $(LPSOLVE) $(OBJ)

$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '\
        o=$(@F) ; b=$${o%.o} ; f=$$b.c ;\
        echo "CC $$f" ;\
        cd $(OBJDIR) ;\
        $(CC) $(CFLAGS) -c -I$(SRCDIR) $(BUNCH_OF_DEFINES) $(SRCDIR)/$$f\
        $(INCLUDE_COLAMD)'
