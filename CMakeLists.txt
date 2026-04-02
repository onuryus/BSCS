cmake_minimum_required(VERSION 3.10)
project(chemsearch)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)

include_directories($ENV{CONDA_PREFIX}/include/rdkit)
link_directories($ENV{CONDA_PREFIX}/lib)

add_executable(build_index src/build_index.cpp)
add_executable(search_index src/search_index.cpp)

target_link_libraries(build_index
    faiss
    RDKitGraphMol
    RDKitSmilesParse
    RDKitRDGeneral
    RDKitFingerprints
    RDKitDataStructs
    OpenMP::OpenMP_CXX
)

target_link_libraries(search_index
    faiss
    RDKitGraphMol
    RDKitSmilesParse
    RDKitRDGeneral
    RDKitFingerprints
    RDKitDataStructs
    OpenMP::OpenMP_CXX
)
