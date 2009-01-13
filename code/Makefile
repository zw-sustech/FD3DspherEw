#******************************************************************************#
#*  Makefile for FD3Dspher program                                            *#
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
MediaMPI  := ON
KernelMPI := ON
#KernelInfoMPI := ON
#KernelAssmMPI := ON

#WITHQS := ON
#DataTypeDouble := ON

SrcSmooth :=ON

#CondFreeCharac := ON
CondFreeTIMG := ON
#CondFreeVHOC := ON
CondFreeVLOW := ON
#CondFreeVEXT := ON

DFLAG_LIST := DEBUG STATIC GETARG VERBOSE \
			  MPITOPO1D USEOMP MPIBARRIER DataTypeDouble \
              WITHQS SrcSmooth  \
			  CondFreeCharac CondFreeTIMG \
              CondFreeVHOC CondFreeVEXT CondFreeVLOW \
              MediaMPI KernelMPI KernelInfoMPI KernelAssmMPI

DFLAGS := $(foreach flag,$(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
DFLAGS := $(strip $(DFLAGS))

#######################################################################
#     list source files and the target names                          #
#######################################################################
skeldirs := OBJ bin checkpoint export input output postproc src
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
EXE_EXPTSNAP   := tool_expt_snap

SRC_EXPTSEISMO := tool_expt_seismo.f90
EXE_EXPTSEISMO := tool_expt_seismo

SRC_KERNEL    := tomo_kernel.f90
EXE_KERNEL    := tomo_kernel
SRC_KERNEL_STA := tomo_kernel_sta.f90
EXE_KERNEL_STA := tomo_kernel_sta
SRC_KERNEL_STA_STRIDE := tomo_kernel_sta_stride.f90
EXE_KERNEL_STA_STRIDE := tomo_kernel_sta_stride
SRC_KERNEL_INFO := tomo_kernel_info.f90
EXE_KERNEL_INFO := tomo_kernel_info
SRC_KERNEL_ASSM := tomo_kernel_assm.f90
EXE_KERNEL_ASSM := tomo_kernel_assm

OBJ_MOD     := $(foreach file,$(SRC_MOD),$(OBJDIR)/$(file:.f90=.o))
OBJ_GRID    := $(foreach file,$(SRC_GRID),$(OBJDIR)/$(file:.f90=.o))
OBJ_MEDIA   := $(foreach file,$(SRC_MEDIA),$(OBJDIR)/$(file:.f90=.o))
OBJ_SOURCE  := $(foreach file,$(SRC_SOURCE),$(OBJDIR)/$(file:.f90=.o))
OBJ_STATION := $(foreach file,$(SRC_STATION),$(OBJDIR)/$(file:.f90=.o))
OBJ_WAVE    := $(foreach file,$(SRC_WAVE),$(OBJDIR)/$(file:.f90=.o))
OBJ_EXPTSNAP   :=  $(foreach file,$(SRC_EXPTSNAP),$(OBJDIR)/$(file:.f90=.o))
OBJ_EXPTSEISMO :=  $(foreach file,$(SRC_EXPTSEISMO),$(OBJDIR)/$(file:.f90=.o))
OBJ_KERNEL     := $(foreach file,$(SRC_KERNEL),$(OBJDIR)/$(file:.f90=.o))
OBJ_KERNEL_STA := $(foreach file,$(SRC_KERNEL_STA),$(OBJDIR)/$(file:.f90=.o))
OBJ_KERNEL_STA_STRIDE := $(foreach file,$(SRC_KERNEL_STA_STRIDE),$(OBJDIR)/$(file:.f90=.o))
OBJ_KERNEL_INFO := $(foreach file,$(SRC_KERNEL_INFO),$(OBJDIR)/$(file:.f90=.o))
OBJ_KERNEL_ASSM := $(foreach file,$(SRC_KERNEL_ASSM),$(OBJDIR)/$(file:.f90=.o))

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

#all: skel preprocess solver postprocess kernel
all: skel preprocess solver kernel

preprocess:  $(BINDIR)/$(EXE_GRID) $(BINDIR)/$(EXE_MEDIA) \
     $(BINDIR)/$(EXE_SOURCE) $(BINDIR)/$(EXE_STATION)
solver:  $(BINDIR)/$(EXE_WAVE)
postprocess: $(BINDIR)/$(EXE_EXPTSNAP) $(BINDIR)/$(EXE_EXPTSEISMO)
kernel: $(BINDIR)/$(EXE_KERNEL) $(BINDIR)/$(EXE_KERNEL_STA) $(BINDIR)/$(EXE_KERNEL_STA_STRIDE) \
     $(BINDIR)/$(EXE_KERNEL_INFO) $(BINDIR)/$(EXE_KERNEL_ASSM)
skel:
	@mkdir -p $(skeldirs)
	@echo "0 0 0 # checkpoint, syncpoint, new nt" > checkpoint.dat

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
$(BINDIR)/$(EXE_KERNEL): $(OBJ_MOD) $(OBJ_KERNEL)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_KERNEL) $(LDFLAGS)
$(BINDIR)/$(EXE_KERNEL_STA): $(OBJ_MOD) $(OBJ_KERNEL_STA)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_KERNEL_STA) $(LDFLAGS)
$(BINDIR)/$(EXE_KERNEL_STA_STRIDE): $(OBJ_MOD) $(OBJ_KERNEL_STA_STRIDE)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_KERNEL_STA_STRIDE) $(LDFLAGS)
$(BINDIR)/$(EXE_KERNEL_INFO): $(OBJ_MOD) $(OBJ_KERNEL_INFO)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_KERNEL_INFO) $(LDFLAGS)
$(BINDIR)/$(EXE_KERNEL_ASSM): $(OBJ_MOD) $(OBJ_KERNEL_ASSM)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_KERNEL_ASSM) $(LDFLAGS)

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
	$(RM) -rf ./output ./checkpoint
	mkdir -p output checkpoint
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

# vim:ft=make:ts=4:sw=4:nu:et:ai:
