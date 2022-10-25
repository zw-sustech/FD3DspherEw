#******************************************************************************#
#*                     Makefile for FD3DspherEw package                       *#
#*                                                                            *#
#*            Author: Wei ZHANG     Email: zhangwei.zw@gmail.com              *#
#*                       Copyright (C) 2008 Wei ZHANG                         *#
#******************************************************************************#

# $Date$
# $Revision$
# $LastChangedBy$

#######################################################################
#                 Compiler and macro flags setting                    #
#######################################################################

#WHEREAMI := nw8440
#WHEREAMI := uri
#WHEREAMI := w520
#WHEREAMI := ustc
#WHEREAMI := generic
WHEREAMI := mars

# Please choose COMPILER in the included Makefile.opt.$(WHEREAMI)

#DEBUG := ON
#STATIC := ON
#GETARG := ON
#VERBOSE :=ON
CheckOverFlow := ON

#MPITOPO1D :=ON
#MPIBARRIER := ON
#USEOMP :=ON
#MPIBuffered := ON

USEPML := ON
#WITHQS := ON
#DataTypeDouble := ON

AnisVTI := ON
#AnisGene := ON

#SrcSurface := ON

#SrcNearestPoint := ON
SrcImpulseFilter := ON
#SrcGaussianSmooth := ON

#CondFreeCharac := ON
CondFreeTIMG := ON
#CondFreeVHOC := ON
CondFreeVLOW := ON
#CondFreeVEXT := ON

DFLAG_LIST := DEBUG STATIC GETARG VERBOSE CheckOverFlow \
			  MPITOPO1D USEOMP MPIBARRIER MPIBuffered DataTypeDouble \
              USEPML WITHQS \
              SrcSurface SrcNearestPoint SrcImpulseFilter SrcGaussianSmooth \
              AnisVTI AnisGene \
			  CondFreeCharac CondFreeTIMG \
              CondFreeVHOC CondFreeVEXT CondFreeVLOW

DFLAGS := $(foreach flag,$(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
DFLAGS := $(strip $(DFLAGS))

#######################################################################
#            directories, source files and target names               #
#######################################################################
skeldirs := OBJ bin
OBJDIR := ./OBJ
BINDIR := ./bin
SRCDIR := ./forward
NRDIR  := ./NR

SRCF90 := mod_constants.F90 mod_math.F90 mod_string.F90 mod_nfseis.F90 \
		  mod_para.F90 mod_mpi.F90  \
		  mod_grid.F90 mod_io.F90 mod_media.F90 mod_macdrp.F90 mod_src.F90 \
		  mod_custom.F90 \
          $(if $(USEPML),mod_abs_npml.F90,mod_abs_exp.F90)
SRCf90 := nrtype.f90 nrutil.f90 bessi0.f90
OBJ_MOD := $(SRCF90:%.F90=%.o) \
          $(SRCf90:%.f90=%.o)
OBJ_MOD := $(addprefix $(OBJDIR)/,$(OBJ_MOD))

EXE_GRID       := seis3d_grid
EXE_METRIC     := seis3d_metric
EXE_MEDIA      := seis3d_media
EXE_MEDIA_MPI  := seis3d_media_mpi
EXE_SOURCE     := seis3d_source
EXE_STATION    := seis3d_station
EXE_WAVE       := seis3d_wave_mpi

EXE_EXPTSNAP   := tool_expt_snap
EXE_EXPTSEISMO := tool_expt_seismo

vpath %.F90 $(SRCDIR)
vpath %.m4  $(SRCDIR)
vpath %.f90 $(NRDIR)


M4       := /usr/bin/m4
M4FLAGS  := $(foreach dir,$(SRCDIR),-I$(dir))

#######################################################################
#                     compiler and option                             #
#######################################################################
include Makefile.opt.$(WHEREAMI)

#######################################################################
#                            target                                   #
#######################################################################
phony_list := skel all preprocess solver postprocess
.PHONY: $(phony_list)

all: skel preprocess solver
#all: skel preprocess solver kernel
#all: skel preprocess solver postprocess kernel

preprocess: $(BINDIR)/$(EXE_GRID)  \
            $(BINDIR)/$(EXE_METRIC)  \
            $(BINDIR)/$(EXE_MEDIA) $(BINDIR)/$(EXE_MEDIA_MPI) \
            $(BINDIR)/$(EXE_SOURCE) $(BINDIR)/$(EXE_STATION)
solver:  $(BINDIR)/$(EXE_WAVE)
#postprocess: $(BINDIR)/$(EXE_EXPTSNAP) $(BINDIR)/$(EXE_EXPTSEISMO)
skel:
	@mkdir -p $(skeldirs)

$(BINDIR)/$(EXE_WAVE): $(OBJ_MOD) $(OBJDIR)/seis3d_wave.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_GRID): $(OBJ_MOD) $(OBJDIR)/seis3d_grid.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_METRIC): $(OBJ_MOD) $(OBJDIR)/seis3d_metric.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_MEDIA): $(OBJ_MOD) $(OBJDIR)/seis3d_media.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_MEDIA_MPI): $(OBJ_MOD) $(OBJDIR)/seis3d_media_mpi.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_SOURCE): $(OBJ_MOD) $(OBJDIR)/seis3d_source.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_STATION): $(OBJ_MOD) $(OBJDIR)/seis3d_station.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_EXPTSNAP): $(OBJ_MOD) $(OBJDIR)/tool_expt_snap.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(BINDIR)/$(EXE_EXPTSEISMO): $(OBJ_MOD) $(OBJDIR)/tool_expt_seismo.o
	$(FC) -o $@ $^ $(LDFLAGS)

RM := rm
cleanexe:
	$(RM) -f $(BINDIR)/*
cleanobj:
	$(RM) -f $(OBJDIR)/*
cleanF90:
	rm -f $(foreach file,$(notdir $(wildcard $(SRCDIR)/*.m4)),$(SRCDIR)/$(file:.m4=.F90))
cleanskel:
	$(RM) -rf $(skeldirs)
	$(RM) -f *.log
#	echo $(foreach file,$(M4DIR)/*,$(FPPDIR)/$(file:.m4=.F90))
#	#$(RM) -f $(foreach file,$(M4DIR)/*,$(FPPDIR)/$(file:.m4=.F90))
cleanall: cleanexe cleanobj cleanF90
distclean: cleanall cleanskel

#######################################################################
#                        subffixes rules                              #
#######################################################################

.SUFFIXES:
.SUFFIXES: .m4 .F90 .f90 .o

$(OBJDIR)/seis3d_media_mpi.o: seis3d_media.F90
	$(FC) $(FFLAGS) $(DFLAGS) -DMediaMPI $< -o $@

%.F90 : %.m4
	$(M4) $(M4FLAGS) $< > $(SRCDIR)/$(@F)

$(OBJDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(NRDIR)/$(<F) -o $@

$(OBJDIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(SRCDIR)/$(<F) -o $@



#######################################################################
#                          include rules                              #
#######################################################################

#-include Makefile.inc.kernel

# vim:ft=make:ts=4:sw=4:nu:et:ai:
