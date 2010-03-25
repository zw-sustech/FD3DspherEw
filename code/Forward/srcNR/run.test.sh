#!/bin/bash

ifort nrtype.f90 nrutil.f90 bessi0.f90 test.f90

./a.out
