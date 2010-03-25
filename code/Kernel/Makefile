#===============================================================================
# Description:
#===============================================================================
#
#                Makefile for SI_kenrel calculation
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
#            macro settings
#===============================================================================

#AnisVTI := ON
#AnisGene := ON
DFLAG_LIST := AnisVTI AnisGene

DFLAGS := $(foreach flag,$(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
DFLAGS := $(strip $(DFLAGS))


#===============================================================================
#            directories, source files and target names
#===============================================================================

#-- forward simulation code module dir
Forward_OBJDIR := ../Forward/OBJ
Forward_SRC_MOD = mod_constants.f90 mod_string.f90 mod_nfseis.f90 \
		  mod_para.f90 mod_mpi.f90  \
		  mod_grid.f90 mod_io.f90 mod_media.f90 
Forward_OBJ_MOD := $(foreach file,$(Forward_SRC_MOD),$(Forward_OBJDIR)/$(file:.f90=.o))

#-- kernel files
skeldirs := OBJ bin src
OBJDIR := ./OBJ
BINDIR := ./bin
SRCDIR := ./src
FPPDIR := ./srcF

#-- source files
SRC_PAIR_FPP := SI_ker_pair.F90
SRC_PAIR := SI_ker_pair.f90
EXE_PAIR := SI_ker_pair
OBJ_PAIR := $(foreach file,$(SRC_PAIR),$(OBJDIR)/$(file:.f90=.o))
SRC_PAIR_MPI := SI_ker_pair_mpi.f90
EXE_PAIR_MPI := SI_ker_pair_mpi
OBJ_PAIR_MPI := $(foreach file,$(SRC_PAIR_MPI),$(OBJDIR)/$(file:.f90=.o))

SRC_STA_FPP := SI_ker_sta.F90
SRC_STA := SI_ker_sta.f90
EXE_STA := SI_ker_sta
OBJ_STA := $(foreach file,$(SRC_STA),$(OBJDIR)/$(file:.f90=.o))
SRC_STA_MPI := SI_ker_sta_mpi.f90
EXE_STA_MPI := SI_ker_sta_mpi
OBJ_STA_MPI := $(foreach file,$(SRC_STA_MPI),$(OBJDIR)/$(file:.f90=.o))

SRC_STA_ANISO_FPP := SI_ker_sta_aniso.F90
SRC_STA_ANISO := SI_ker_sta_aniso.f90
EXE_STA_ANISO := SI_ker_sta_aniso
OBJ_STA_ANISO := $(foreach file,$(SRC_STA_ANISO),$(OBJDIR)/$(file:.f90=.o))
SRC_STA_ANISO_MPI := SI_ker_sta_aniso_mpi.f90
EXE_STA_ANISO_MPI := SI_ker_sta_aniso_mpi
OBJ_STA_ANISO_MPI := $(foreach file,$(SRC_STA_ANISO_MPI),$(OBJDIR)/$(file:.f90=.o))

SRC_KER2INFO_FPP := SI_ker2info.F90
SRC_KER2INFO := SI_ker2info.f90
EXE_KER2INFO := SI_ker2info
OBJ_KER2INFO := $(foreach file,$(SRC_KER2INFO),$(OBJDIR)/$(file:.f90=.o))
SRC_KER2INFO_MPI := SI_ker2info_mpi.f90
EXE_KER2INFO_MPI := SI_ker2info_mpi
OBJ_KER2INFO_MPI := $(foreach file,$(SRC_KER2INFO_MPI),$(OBJDIR)/$(file:.f90=.o))

SRC_KER2ASSM_FPP := SI_ker2assm.F90
SRC_KER2ASSM := SI_ker2assm.f90
EXE_KER2ASSM := SI_ker2assm
OBJ_KER2ASSM := $(foreach file,$(SRC_KER2ASSM),$(OBJDIR)/$(file:.f90=.o))
SRC_KER2ASSM_MPI := SI_ker2assm_mpi.f90
EXE_KER2ASSM_MPI := SI_ker2assm_mpi
OBJ_KER2ASSM_MPI := $(foreach file,$(SRC_KER2ASSM_MPI),$(OBJDIR)/$(file:.f90=.o))

vpath %.F90 $(FPPDIR)

#===============================================================================
#                 Compiler and macro flags setting
#===============================================================================
FC     :=  /opt/mpich/p4-intel/bin/mpif90
NETCDF :=  /opt/netcdf/intel-10.1

LDFLAGS := -L$(NETCDF)/lib -lnetcdf \
           -Bstatic -static
FFLAGS := -c -O3 -w95 -warn all -zero \
          -module $(OBJDIR) -module $(Forward_OBJDIR) \
          -I$(FPPDIR) \
          -I$(NETCDF)/include 

FPP      := /usr/bin/cpp
FPPFLAGS := -P -traditional $(foreach dir,$(FPPDIR),-I$(dir)) $(DFLAGS) 

#===============================================================================
#            target
#===============================================================================
phony_list := kernel
.PHONY: $(phony_list)

all: skel kernel

skel:
	@mkdir -p $(skeldirs)

#kernel: $(BINDIR)/$(EXE_STA_ANISO)  $(BINDIR)/$(EXE_STA_ANISO_MPI)
#kernel: $(BINDIR)/$(EXE_KER2ASSM) $(BINDIR)/$(EXE_KER2ASSM_MPI)

kernel: $(BINDIR)/$(EXE_STA)  $(BINDIR)/$(EXE_STA_MPI) \
        $(BINDIR)/$(EXE_STA_ANISO)  $(BINDIR)/$(EXE_STA_ANISO_MPI) \
        $(BINDIR)/$(EXE_KER2ASSM) $(BINDIR)/$(EXE_KER2ASSM_MPI)

#kernel: $(BINDIR)/$(EXE_PAIR) $(BINDIR)/$(EXE_PAIR_MPI) \
#         $(BINDIR)/$(EXE_STA)  $(BINDIR)/$(EXE_STA_MPI) \
#        $(BINDIR)/$(EXE_KER2INFO) $(BINDIR)/$(EXE_KER2INFO_MPI) \
#        $(BINDIR)/$(EXE_KER2ASSM) $(BINDIR)/$(EXE_KER2ASSM_MPI)

$(BINDIR)/$(EXE_PAIR): $(Forward_OBJ_MOD) $(OBJ_PAIR)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_PAIR) $(LDFLAGS)
$(BINDIR)/$(EXE_PAIR_MPI): $(Forward_OBJ_MOD) $(OBJ_PAIR_MPI)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_PAIR_MPI) $(LDFLAGS)

$(BINDIR)/$(EXE_STA): $(Forward_OBJ_MOD) $(OBJ_STA)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_STA) $(LDFLAGS)
$(BINDIR)/$(EXE_STA_MPI): $(Forward_OBJ_MOD) $(OBJ_STA_MPI)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_STA_MPI) $(LDFLAGS)

$(BINDIR)/$(EXE_STA_ANISO): $(Forward_OBJ_MOD) $(OBJ_STA_ANISO)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_STA_ANISO) $(LDFLAGS)
$(BINDIR)/$(EXE_STA_ANISO_MPI): $(Forward_OBJ_MOD) $(OBJ_STA_ANISO_MPI)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_STA_ANISO_MPI) $(LDFLAGS)

$(BINDIR)/$(EXE_KER2INFO): $(Forward_OBJ_MOD) $(OBJ_KER2INFO)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_KER2INFO) $(LDFLAGS)
$(BINDIR)/$(EXE_KER2INFO_MPI): $(Forward_OBJ_MOD) $(OBJ_KER2INFO_MPI)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_KER2INFO_MPI) $(LDFLAGS)

$(BINDIR)/$(EXE_KER2ASSM): $(Forward_OBJ_MOD) $(OBJ_KER2ASSM)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_KER2ASSM) $(LDFLAGS)
$(BINDIR)/$(EXE_KER2ASSM_MPI): $(Forward_OBJ_MOD) $(OBJ_KER2ASSM_MPI)
	$(FC) -o $@ $(Forward_OBJ_MOD) $(OBJ_KER2ASSM_MPI) $(LDFLAGS)

RM := rm
cleanexe:
	$(RM) -f $(BINDIR)/*
cleanobj:
	$(RM) -f $(OBJDIR)/*
cleanf90:
	$(RM) -f $(SRCDIR)/*
cleanskel:
	$(RM) -rf $(skeldirs)
	$(RM) -f *.log
cleanall: cleanexe cleanobj cleanf90 
distclean: cleanall cleanskel

#######################################################################
#                        subffixes rules                              #
#######################################################################
$(SRC_PAIR_MPI): $(SRC_PAIR_FPP)
	$(FPP) $(FPPINCS) $(FPPFLAGS) -DKernelMPI $< > $(SRCDIR)/$(@F)
$(SRC_STA_MPI): $(SRC_STA_FPP)
	$(FPP) $(FPPINCS) $(FPPFLAGS) -DKernelMPI $< > $(SRCDIR)/$(@F)
$(SRC_STA_ANISO_MPI): $(SRC_STA_ANISO_FPP)
	$(FPP) $(FPPINCS) $(FPPFLAGS) -DKernelMPI $< > $(SRCDIR)/$(@F)
$(SRC_KER2INFO_MPI): $(SRC_KER2INFO_FPP)
	$(FPP) $(FPPINCS) $(FPPFLAGS) -DKernelInfoMPI $< > $(SRCDIR)/$(@F)
$(SRC_KER2ASSM_MPI): $(SRC_KER2ASSM_FPP)
	$(FPP) $(FPPINCS) $(FPPFLAGS) -DKernelAssmMPI $< > $(SRCDIR)/$(@F)

%.f90 : %.F90
	$(FPP) $(FPPFLAGS) $< > $(SRCDIR)/$(@F)

$(OBJDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(SRCDIR)/$(<F) -o $@

# vim:ft=make:ts=4:sw=4:nu:et:ai:
