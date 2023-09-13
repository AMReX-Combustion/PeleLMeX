#!/usr/bin/env bash
#
# Copyright 2020 Axel Huebl
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential     \
    ca-certificates     \
    cmake               \
    gnupg               \
    libopenmpi-dev      \
    openmpi-bin         \
    pkg-config          \
    wget

curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb
sudo dpkg -i cuda-keyring_1.0-1_all.deb
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-11-2 \
    cuda-compiler-11-2           \
    cuda-cupti-dev-11-2          \
    cuda-minimal-build-11-2      \
    cuda-nvml-dev-11-2           \
    cuda-nvtx-11-2               \
    libcurand-dev-11-2           \
    libcusparse-dev-11-2         \
    libcusolver-dev-11-2         \
    libcublas-dev-11-2

sudo ln -s cuda-11.2 /usr/local/cuda
