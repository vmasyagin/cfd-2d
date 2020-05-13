cmake_minimum_required(VERSION 3.8)
project(cfd_2d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CONTRIB_DIR "${CMAKE_SOURCE_DIR}/contrib")

set(CONFIGS_DIR "${CMAKE_SOURCE_DIR}/configs")

configure_file(config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h HEADER_FILE_ONLY)

# необходимо прописать путь к библиотеке MPI
set(MPI_INCLUDE "${CONTRIB_DIR}/MPICH2/include")
set(MPI_LIB "${CONTRIB_DIR}/MPICH2/lib")

#
#  http://acts.nersc.gov/hypre/
#
set(HYPRE_INCLUDE "${CONTRIB_DIR}/HYPRE/include")
set(HYPRE_LIB "${CONTRIB_DIR}/HYPRE/lib")


#
#  http://glaros.dtc.umn.edu/gkhome/metis/metis/download
#
set(METIS_INCLUDE "${CONTRIB_DIR}/METIS/include")
set(METIS_LIB "${CONTRIB_DIR}/METIS/lib")

#
#  http://glaros.dtc.umn.edu/gkhome/metis/metis/download
#
set(AMGX_INCLUDE "${CONTRIB_DIR}/AMGX/include")
set(AMGX_LIB "${CONTRIB_DIR}/AMGX/lib")

include_directories(
		${MPI_INCLUDE}
		${HYPRE_INCLUDE}
		${METIS_INCLUDE}
		${AMGX_INCLUDE}
		src
		src/methods
		src/methods/solvers
		src/mesh
		src/tinyxml)
link_directories(
		${MPI_LIB}
		${HYPRE_LIB}
		${METIS_LIB}
		${AMGX_LIB}
)

set(METHODS_SRC
		src/methods/decomp.cpp
		src/methods/solver.cpp
		src/methods/fvm_tvd.cpp
		src/methods/fvm_tvd_implicit.cpp
		src/methods/fem_rkdg.cpp
		src/methods/fem_dg_implicit.cpp
		src/methods/LimiterDG.cpp
		src/methods/LimiterDGCockburn.cpp
)

set(SOLVERS_SRC
		src/methods/solvers/CSR.cpp
		src/methods/solvers/MatrixSolver.cpp
		src/methods/solvers/SolverHypre.cpp
		src/methods/solvers/SolverHypreBoomerAmg.cpp
		src/methods/solvers/SolverHypreCustomSeidel.cpp
		src/methods/solvers/SolverHypreFlexGmres.cpp
		src/methods/solvers/SolverHypreFlexGmresPrecAMG.cpp
		src/methods/solvers/SolverHypreGmres.cpp
		src/methods/solvers/SolverHyprePcg.cpp
		src/methods/solvers/SolverZeidel.cpp
		src/methods/solvers/SolverAMGX.cpp
		src/methods/solvers/SolverAMGXImpl.cpp
        src/methods/solvers/Initializator.h src/methods/solvers/Initializator.cpp)

set(MESH_SRC
		src/mesh/grid.cpp
		src/mesh/MeshReader.cpp
		src/mesh/MeshReaderBerkleyTriangle.cpp
		src/mesh/MeshReaderSalomeUnv.cpp
		src/mesh/MeshReaderGmsh.cpp
)

set(TINYXML_SRC
		src/tinyxml/tinystr.cpp
		src/tinyxml/tinyxml.cpp
		src/tinyxml/tinyxmlerror.cpp
		src/tinyxml/tinyxmlparser.cpp
)

add_executable(	cfd-2d
		src/main.cpp
		src/global.cpp
		src/bnd_cond.cpp
		${SOLVERS_SRC}
		${METHODS_SRC}
		${MESH_SRC}
		${TINYXML_SRC}
)

target_link_libraries(cfd-2d mpi amgxsh HYPRE metis)
#install(TARGETS cfd-2d RUNTIME DESTINATION bin)
