cmake_minimum_required( VERSION 2.6 )

set( CMAKE_BUILD_TYPE Release )
set( CMAKE_CXX_STANDARD 11 )

set( PATH_TO_PROJECT /Users/wuelker/Documents/ProteinDocking )
set( PATH_TO_LIBRARIES /usr/local/lib )
set( PATH_TO_INCLUDE /usr/local/include )

project( ProteinDocking )

add_definitions( -DFIGTREE_NO_ANN )

add_executable( ProteinDocking main.cpp ./include/SpharmonicKit27/cospmls.c ./include/SpharmonicKit27/primitive.c ./include/SpharmonicKit27/newFCT.c ./include/SpharmonicKit27/OURperms.c ./include/SpharmonicKit27/OURmods.c ./include/SpharmonicKit27/FST_semi_memo.c ./include/SpharmonicKit27/fft_grids.c ./include/SpharmonicKit27/naive_synthesis.c ./include/SpharmonicKit27/primitive_FST.c ./include/SpharmonicKit27/seminaive.c ./include/SpharmonicKit27/csecond.c ./include/SpharmonicKit27/oddweights.c ./include/SpharmonicKit27/weights.c ./include/SpharmonicKit27/FFTcode.c ./include/SpharmonicKit27/indextables.c ./include/SpharmonicKit27/permroots.c ./include/figtree-master/src/figtree.cpp ./include/figtree-master/src/KCenterClustering.cpp ./include/fgt.cpp ./include/fsft.cpp ./include/fsglft.cpp )

link_directories( ${PATH_TO_LIBRARIES} )
target_link_libraries( ProteinDocking nfft3 )

include_directories( ${PATH_TO_PROJECT}/include )
include_directories( ${PATH_TO_PROJECT}/include/figtree-master/include )
include_directories( ${PATH_TO_INCLUDE} )
