# billion-scale-chemical-search

# Billion-Scale Chemical Similarity Search

A high-performance chemical search engine built with C++, RDKit and FAISS for million to billion scale molecular libraries.

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
```

## 🚀 Usage

### 1. Prepare data

Input file must be a SMILES file named:

```bash
data/tam.smi
./build_index
./search_index
SMILES (exit ile cik): CC1=CC=CC=C1
```

<table>
  <thead>
    <tr>
      <th>Parameter</th>
      <th>Description</th>
      <th>Effect on Accuracy</th>
      <th>Effect on Speed</th>
      <th>Default</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>--metric</code></td>
      <td>Similarity metric used for reranking results</td>
      <td>High (defines final ranking)</td>
      <td>Low impact</td>
      <td><code>tanimoto</code></td>
    </tr>
    <tr>
      <td><code>--k</code></td>
      <td>Number of candidates retrieved from FAISS before reranking</td>
      <td>Higher → better recall</td>
      <td>Higher → slower reranking</td>
      <td><code>5000</code></td>
    </tr>
    <tr>
      <td><code>--nprobe</code></td>
      <td>Number of clusters searched in FAISS (IVF search scope)</td>
      <td>Higher → better accuracy</td>
      <td>Higher → slower search</td>
      <td><code>32</code></td>
    </tr>
    <tr>
      <td><code>--max_codes</code></td>
      <td>Maximum number of codes scanned during search (limits computation)</td>
      <td>Lower → less accurate</td>
      <td>Lower → faster</td>
      <td><code>0</code> (disabled)</td>
    </tr>
  </tbody>
</table>
