#!/usr/bin/env bash

set -eu

datadir=/scratch/snx3000/juckerj/WG_1_tasks/ecrad/models
exedir=/scratch/snx3000/juckerj/WG_1_tasks/ecrad/infero/ecrad


eng_type=tf_c

model_path=$datadir/serialized_model_step_34000.tf

# run inference
$exedir/real_model.exe \
  $model_path \
  $eng_type \


