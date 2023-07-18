#!/bin/bash

# a script based on https://gist.github.com/nicholst/009b22b9a6f94d93ed50fb27ea6f0f82
# to visualize the results of a GLM analysis in FSLeyes

EffectImg=/Users/alexandresayal/GitHub/vpinhibitioncibit/data/expC-glm-sub01/spmT_0001_thres_fwe0.01_k100.img
StatImg=/Users/alexandresayal/GitHub/vpinhibitioncibit/data/expC-glm-sub01/spmT_0001.nii

# Value to set outline of signficant voxels
PosStatThresh=5.50
NegStatThresh=-999  # Set to large negative value to have no -ev threshold

# Color of outline
StatThreshColor="0.2 0.2 0.2"  # I find this dark gray 20% value work well

# Effect will be opaque when stat value exceeds this value
StatMaxModVal=15

fsleyes \
  --showColourBar --colourBarLocation right --colourBarLabelSide top-left \
  -std \
  $EffectImg \
   --cmap red-yellow --negativeCmap blue-lightblue --useNegativeCmap \
   --modulateAlpha --modulateImage $StatImg --modulateRange 0.0 $StatMaxModVal \
  $StatImg --overlayType mask --outline \
   --maskColour $StatThreshColor --threshold $NegStatThresh $PosStatThresh