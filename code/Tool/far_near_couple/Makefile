#===============================================================================
# Description:
#===============================================================================
#
#                Makefile for tool_coupling_preprocess.F90
# 
#            Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
#                       Copyright (C) 2010 Wei ZHANG
#
#-------------------------------------------------------------------------------
# Time stamp, log, version:
#-------------------------------------------------------------------------------
#
# $Date$
# $Revision$
# $LastChangedBy$
#

#===============================================================================
#                 Compiler and macro flags setting
#===============================================================================

FC     :=  /opt/mpich/p4-intel/bin/mpif90
NETCDF :=  /opt/netcdf/intel-10.1

#===============================================================================
#            directories, source files and target names
#===============================================================================

#-- main code module dir
OBJDIR := ../../OBJ
FPPDIR := ../../srcF
SRCDIR := ../../src

#-- activate DEBUG macro to print debug information in the code
#DFLAGS := -DDEBUG

FPP      := /usr/bin/cpp
FPPFLAGS := -P -traditional $(foreach dir,$(FPPDIR),-I$(dir)) $(DFLAGS) 

M4       := /usr/bin/m4
M4FLAGS  := $(foreach dir,$(M4DIR),-I$(dir))

LDFLAGS := -L$(NETCDF)/lib -lnetcdf \
           -Bstatic -static
FFLAGS := -c -O3 -w95 -warn all -zero -module $(OBJDIR) -I$(NETCDF)/include 

#vpath %.F90 $(FPPDIR)
#vpath %.m4  $(M4DIR)

#===============================================================================
#            target
#===============================================================================
SRC := tool_coupling_preprocess.f90
EXE := tool_coupling_preprocess
OBJ := $(foreach file,$(SRC),$(file:.f90=.o))

SRC_MOD = mod_constants.f90 mod_string.f90 mod_nfseis.f90 \
		  mod_para.f90 mod_mpi.f90  \
		  mod_grid.f90 mod_io.f90 mod_media.f90 mod_macdrp.f90 mod_src.f90
OBJ_MOD := $(foreach file,$(SRC_MOD),$(OBJDIR)/$(file:.f90=.o))

$(EXE): $(OBJ)
	$(FC) -o $@ $(OBJ) $(OBJ_MOD) $(LDFLAGS)

RM := rm
clean:
	$(RM) -f $(EXE) *.o *.mod

#######################################################################
#                        subffixes rules                              #
#######################################################################

.SUFFIXES:
.SUFFIXES: .m4 .F90 .f90 .o

%.f90 : %.m4
	$(M4) $(M4FLAGS) $< > $(FPPDIR)/$(@F:.f90=.F90)
	$(FPP) $(FPPFLAGS) $(FPPDIR)/$(@F:.f90=.F90) > $(SRCDIR)/$(@F)

%.f90 : %.F90
	$(FPP) $(FPPFLAGS) $< > $(@F)

%.o : %.f90
	$(FC) $(FFLAGS) $(<F) -o $@

# vim:ft=make:ts=4:sw=4:nu:et:ai:
