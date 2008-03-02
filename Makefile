#******************************************************************************#
#*  Makefile for FD3DTOPO-macdrp-spher program                                *#
#*                                                                            *#
#*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *#
#*  Copyright (C) Wei ZHANG, 2007. All Rights Reserved.                       *#
#******************************************************************************#

# $Date$
# $Revision$
# $LastChangedBy$

#######################################################################
#                 Compiler and macro flags setting                    #
#######################################################################

#WHEREAMI := nw8440
#WHEREAMI := x60s
#WHEREAMI := pku
#WHEREAMI := ustc
WHEREAMI := uri
#WHEREAMI := generic

# Please choose COMPILER in the included Makefile.$(WHEREAMI)

#DEBUG := ON
STATIC := ON
#GETARG := ON
#VERBOSE :=ON

#MPITOPO1D :=ON
#MPIBARRIER := ON
#USEOMP :=ON

#WITHQS := ON
#WithoutVHOC := ON
#DataTypeDouble := ON

#SrcSmooth :=ON
#SrcTensorMomentum := ON
SrcTensorHook := ON
SrcForce := ON

#CondFreeCharac := ON
CondFreeTIMG := ON

DFLAG_LIST := DEBUG STATIC GETARG VERBOSE \
			  MPITOPO1D USEOMP MPIBARRIER \
              WITHQS WithoutVHOC DataTypeDouble \
			  SrcSmooth SrcTensorMomentum SrcTensorHook SrcForce \
			  CondFreeCharac CondFreeTIMG

DFLAGS := $(foreach flag,$(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
DFLAGS := $(strip $(DFLAGS))

#######################################################################
#     list source files and the target names                          #
#######################################################################
skeldirs := OBJ bin export input output postproc rest src
OBJDIR := ./OBJ
BINDIR := ./bin
SRCDIR := ./src
FPPDIR := ./srcF

SRC_MOD = mod_constants.f90 mod_math.f90 mod_string.f90 mod_nfseis.f90 \
		  mod_para.f90 mod_mpi.f90  \
		  mod_grid.f90 mod_io.f90 mod_media.f90 mod_macdrp.f90 mod_src.f90 \
		  mod_custom.f90
SRC_WAVE := $(if $(USEPML),mod_abs_npml.f90,mod_abs_exp.f90) \
          seis3d_wave.f90
SRC_GRID    := seis3d_grid.f90
SRC_MEDIA   := seis3d_media.f90
SRC_SOURCE  := seis3d_source.f90
SRC_STATION := seis3d_station.f90
EXE_GRID    := seis3d_grid
EXE_MEDIA   := seis3d_media
EXE_SOURCE  := seis3d_source
EXE_STATION := seis3d_station
EXE_WAVE    := seis3d_wave

SRC_EXPTSNAP   := tool_expt_snap.f90
SRC_EXPTSEISMO := tool_expt_seismo.f90
EXE_EXPTSNAP   := tool_expt_snap
EXE_EXPTSEISMO := tool_expt_seismo


OBJ_MOD     := $(foreach file,$(SRC_MOD),$(OBJDIR)/$(file:.f90=.o))
OBJ_GRID    := $(foreach file,$(SRC_GRID),$(OBJDIR)/$(file:.f90=.o))
OBJ_MEDIA   := $(foreach file,$(SRC_MEDIA),$(OBJDIR)/$(file:.f90=.o))
OBJ_SOURCE  := $(foreach file,$(SRC_SOURCE),$(OBJDIR)/$(file:.f90=.o))
OBJ_STATION := $(foreach file,$(SRC_STATION),$(OBJDIR)/$(file:.f90=.o))
OBJ_WAVE    := $(foreach file,$(SRC_WAVE),$(OBJDIR)/$(file:.f90=.o))
OBJ_EXPTSNAP   :=  $(foreach file,$(SRC_EXPTSNAP),$(OBJDIR)/$(file:.f90=.o))
OBJ_EXPTSEISMO :=  $(foreach file,$(SRC_EXPTSEISMO),$(OBJDIR)/$(file:.f90=.o))

vpath %.F90 $(FPPDIR)

FPP := /usr/bin/cpp
FPPFLAGS := -P -traditional $(foreach dir,$(FPPDIR),-I$(dir)) $(DFLAGS) 

#######################################################################
#                     compiler and option                             #
#######################################################################
include Makefile.$(WHEREAMI)

#######################################################################
#                        the target                                   #
#######################################################################
.PHONY: skel all preprocess solver postprocess kernel

all: skel preprocess solver postprocess kernel

-include Makefile.kernel

preprocess:  $(BINDIR)/$(EXE_GRID) $(BINDIR)/$(EXE_GRID_MPI)  \
     $(BINDIR)/$(EXE_MEDIA) $(BINDIR)/$(EXE_MEDIA_MPI) \
     $(BINDIR)/$(EXE_SOURCE) $(BINDIR)/$(EXE_STATION)
solver:  $(BINDIR)/$(EXE_WAVE)
postprocess: $(BINDIR)/$(EXE_EXPTSNAP) $(BINDIR)/$(EXE_EXPTSEISMO)
kernel: $(BINDIR)/$(EXE_KERNEL)
skel:
	@mkdir -p $(skeldirs)
	@echo 0 > rest_point.dat

$(BINDIR)/$(EXE_WAVE): $(OBJ_MOD)  $(OBJ_WAVE)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_WAVE) $(LDFLAGS)
$(BINDIR)/$(EXE_GRID): $(OBJ_MOD) $(OBJ_GRID)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_GRID) $(LDFLAGS)
$(BINDIR)/$(EXE_MEDIA): $(OBJ_MOD) $(OBJ_MEDIA)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_MEDIA) $(LDFLAGS)
$(BINDIR)/$(EXE_SOURCE): $(OBJ_MOD) $(OBJ_SOURCE)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_SOURCE) $(LDFLAGS)
$(BINDIR)/$(EXE_STATION): $(OBJ_MOD) $(OBJ_STATION)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_STATION) $(LDFLAGS)
$(BINDIR)/$(EXE_EXPTSNAP): $(OBJ_MOD) $(OBJ_EXPTSNAP)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_EXPTSNAP) $(LDFLAGS)
$(BINDIR)/$(EXE_EXPTSEISMO): $(OBJ_MOD) $(OBJ_EXPTSEISMO)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_EXPTSEISMO) $(LDFLAGS)

RM := rm
cleanexe:
	$(RM) -f $(BINDIR)/*
cleanobj:
	$(RM) -f $(OBJDIR)/*
cleanf90:
	$(RM) -f $(SRCDIR)/*
cleaninput:
	$(RM) -rf ./input/
	mkdir -p input
cleanoutput:
	$(RM) -rf ./output ./rest/
	mkdir -p output rest
cleanall: cleanexe cleanobj cleanf90
distclean: cleanexe cleanobj cleanf90 cleaninput cleanoutput

#######################################################################
#                        subffixes rules                              #
#######################################################################

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

%.f90 : %.F90
	$(FPP) $(FPPFLAGS) $< > $(SRCDIR)/$(@F)

$(OBJDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(SRCDIR)/$(<F) -o $@
