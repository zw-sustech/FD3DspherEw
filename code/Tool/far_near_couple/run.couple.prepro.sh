#!/bin/bash

EXE=./tool_coupling_preprocess

#ntwin_Far_field=6000
ntwin_Far_field=10

is_couple_at_x1='.true.'
index_at_x1="27 32 27 132 27 130"
fnm_x1="./wave_far_x1.bin"

#is_couple_at_x2='.false.'
is_couple_at_x2='.true.'
#index_at_x2="0 0 0 0 0 0"
index_at_x2="127 132 27 132 27 130"
fnm_x2="./wave_far_x2.bin"

#is_couple_at_y1='.false.'
is_couple_at_y1='.true.'
index_at_y1="27 132 27 32 27 130"
fnm_y1="./wave_far_y1.bin"

#is_couple_at_y2='.false.'
is_couple_at_y2='.true.'
index_at_y2="27 132 127 132 27 130"
fnm_y2="./wave_far_y2.bin"

#is_couple_at_z1='.false.'
is_couple_at_z1='.true.'
index_at_z1="27 132 27 132 27 32"
fnm_z1="./wave_far_z1.bin"

is_couple_at_z2='.false.'
index_at_z2="0 0 0 0 0 0"
fnm_z2="none"

output_dir="./input"

$EXE << -EOF-
${ntwin_Far_field}
$is_couple_at_x1
$index_at_x1
$fnm_x1
$is_couple_at_x2
$index_at_x2
$fnm_x2
$is_couple_at_y1
$index_at_y1
$fnm_y1
$is_couple_at_y2
$index_at_y2
$fnm_y2
$is_couple_at_z1
$index_at_z1
$fnm_z1
$is_couple_at_z2
$index_at_z2
$fnm_z2
$output_dir
-EOF-
