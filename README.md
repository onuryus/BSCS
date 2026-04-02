# billion-scale-chemical-search

# Billion-Scale Chemical Similarity Search

A high-performance chemical search engine built with RDKit and FAISS for million- to billion-scale molecular libraries.

It uses Morgan fingerprints, FAISS IVFPQ indexing, and exact reranking with multiple similarity metrics.


## Features

- Morgan fingerprint generation with RDKit
- FAISS IVFPQ approximate nearest neighbor search
- Streaming index construction for large libraries
- Offset-based SMILES retrieval
- Exact reranking with:
  - Tanimoto
  - Dice
  - Tversky
  - Cosine
  - Kulczynski
- Search-time tuning with:
  - `--metric`
  - `--k`
  - `--nprobe`
  - `--max_codes`
 
## How it works

text
SMILES
  -> RDKit molecule
  -> Morgan fingerprint (2048-bit)
  -> float vector conversion
  -> FAISS IVFPQ indexing
  -> approximate candidate retrieval
  -> exact similarity reranking



## ⚙️ Installation

### 1. Create environment

```bash
mamba create -n chem_cpp \
  -c conda-forge \
  python=3.10 \
  cmake make gcc gxx \
  boost-cpp eigen \
  faiss-cpu \
  -y

conda activate chem_cpp



git clone https://github.com/rdkit/rdkit.git
cd rdkit
mkdir build && cd build

cmake \
  -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
  -DRDK_INSTALL_INTREE=OFF \
  -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
  -DRDK_BUILD_CAIRO_SUPPORT=OFF \
  -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
  ..

make -j$(nproc)
make install




git clone https://github.com/YOUR_USERNAME/chem-faiss-search.git
cd chem-faiss-search

rm -rf build
mkdir build && cd build

cmake ..
make -j$(nproc)
