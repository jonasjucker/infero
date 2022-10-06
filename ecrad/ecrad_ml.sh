#!/usr/bin/env bash

set -eu

datadir=$1
eng_type=tf_c

model=$datadir/patch_shrNhbrWeightslvl0_all_heights_BN_mse_datav2IG_SSD_MT_featNormv2_fit_seq_MergedLwupLwdnSwupSwdn/serialized_model_steps_182400_238640_80560_188480

# run inference
./multiple_fields_model.exe \
  $model \
  $eng_type \


