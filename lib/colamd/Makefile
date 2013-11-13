#--------------------------------------------------#
# Makefile to create a static library for LPSOLVE. #
#                                                  #
# PLR 11.2013                                      #
#--------------------------------------------------#

# Basic stuff
SRCDIR     = $(PWD)
OBJDIR     = $(SRCDIR)/obj
LIBDIR     = $(SRCDIR)/..
COLAMD_AR  = libcolamd.a
COLAMD     = $(LIBDIR)/$(COLAMD_AR)
.SILENT:

# Defaults for system-dependent variables
AR     = ar crs
CC     = gcc
CFLAGS = -O3

# Main target
.PHONY: default
default: $(COLAMD)

# Target for deleting object files
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(COLAMD)

# Object-file list
OBJ = \
$(OBJDIR)/colamd.o
# End object-file list

# How to compile things
$(COLAMD): $(OBJ)
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        echo "AR COLAMD" ;\
        rm -f $(COLAMD) ;\
        $(AR) $(COLAMD) $(OBJ)

$(OBJ):
	[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '\
        o=$(@F) ; b=$${o%.o} ; f=$$b.c ;\
        echo "CC $$f" ;\
        cd $(OBJDIR) ;\
        $(CC) $(CFLAGS) -c -I$(SRCDIR) $(SRCDIR)/$$f'
