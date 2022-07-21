#!/bin/bash

# # initialise module environment if it is not
# if [[ ! $(command -v module > /dev/null 2>&1) ]]; then
#   . /usr/local/apps/module/init/bash
# fi
# module load cmake

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# ======= Basic config =========
INFERO_VERBOSE_COMPILATION=${INFERO_VERBOSE_COMPILATION:-"0"}
ROOT_DIR=${INFERO_HOME:-$HOME}

ROOT_SRC_DIR=${ROOT_DIR}/local
ROOT_BUILD_DIR=${ROOT_DIR}/builds
ROOT_INSTALL_DIR=${ROOT_DIR}/installs

# MPI support
WITH_MPI=OFF

# Fortran API support
WITH_FCKIT=ON

# ML backends
WITH_ONNX_RUNTIME=OFF
WITH_TFC_RUNTIME=OFF
WITH_TFLITE_RUNTIME=OFF
WITH_TRT=OFF

# Tests
ENABLE_TESTS=ON

# build #procs
BUILD_NPROCS=8
# ===============================


# ======== other configs ========

# ECBUILD
ECBUILD_BRANCH=develop
ECBUILD_SRC_DIR=${ROOT_SRC_DIR}/ecbuild
ECBUILD_BUILD_DIR=${ECBUILD_SRC_DIR}
ECBUILD_BUILD_EXE=${ECBUILD_BUILD_DIR}/bin/ecbuild

# ECKIT
ECKIT_BRANCH=1.20.0
ECKIT_SRC_DIR=${ROOT_SRC_DIR}/eckit
ECKIT_BUILD_DIR=${ROOT_BUILD_DIR}/eckit

# FCKIT
FCKIT_BRANCH=develop
FCKIT_SRC_DIR=${ROOT_SRC_DIR}/fckit
FCKIT_BUILD_DIR=${ROOT_BUILD_DIR}/fckit

# ONNX runtime
arch=$(uname -m)
if [[ "${OSTYPE}" == "linux"* ]] && [[ "${arch}" == "x86_64" ]]; then
  ONNX_VERSION=1.10.0
  ONNX_TARFILE=onnxruntime-linux-x64-${ONNX_VERSION}.tgz
elif [[ "${OSTYPE}" == "darwin"* ]] && [[ "${arch}" == "arm64" ]]; then
  ONNX_VERSION=1.11.1
  ONNX_TARFILE=onnxruntime-osx-arm64-${ONNX_VERSION}.tgz
fi
ONNX_SOURCE_DIR=${ROOT_SRC_DIR}/onnxruntime
ONNX_BUILD_DIR=${ROOT_SRC_DIR}/onnxruntime
ONNX_URL=https://github.com/microsoft/onnxruntime/releases/download
    
# TF_C
TFC_VERSION=2.6.0
TFC_SOURCE_DIR=${ROOT_SRC_DIR}/tensorflow_c
TFC_BUILD_DIR=${ROOT_SRC_DIR}/tensorflow_c
TFC_URL=https://storage.googleapis.com/tensorflow/libtensorflow
TFC_TARFILE=libtensorflow-cpu-linux-x86_64-${TFC_VERSION}.tar.gz

# TFLITE
TFLITE_SOURCE_DIR=${ROOT_SRC_DIR}/tensorflow
TFLITE_BUILD_DIR=${ROOT_BUILD_DIR}/tflite

# TENSORRT
# NB TensorRT must be downloaded separately..
TRT_SOURCE_DIR=${ROOT_SRC_DIR}/TensorRT-8.0.0.3
TRT_BUILD_DIR=${ROOT_SRC_DIR}/TensorRT-8.0.0.3

# INFERO
INFERO_SRC_DIR=$(dirname ${SCRIPT_DIR})
INFERO_BUILD_DIR=${ROOT_BUILD_DIR}/infero
# ===============================

