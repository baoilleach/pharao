cmake_minimum_required (VERSION 2.6) 

# Project
project(pharao)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
find_package(OpenBabel2 REQUIRED)

# Version
set (${PROJECT_NAME}_VERSION_MAJOR 3)
set (${PROJECT_NAME}_VERSION_MINOR 0)
set (${PROJECT_NAME}_VERSION_PATCH 3)
configure_file ("${PROJECT_SOURCE_DIR}/pharaoConfig.h.in" "${PROJECT_SOURCE_DIR}/include/pharaoConfig.h")

# Release type
set (CMAKE_BUILD_TYPE Release)

# Sieve code
set (PharaoDir ${PROJECT_SOURCE_DIR}/src)
set (sources     ${PharaoDir}/main.cpp
                 ${PharaoDir}/addBest.cpp
                 ${PharaoDir}/alignment.cpp
                 ${PharaoDir}/aromFuncCalc.cpp
                 ${PharaoDir}/calcPharm.cpp
                 ${PharaoDir}/chargeFuncCalc.cpp
                 ${PharaoDir}/compScore.cpp
                 ${PharaoDir}/coordinate.cpp
                 ${PharaoDir}/functionMapping.cpp
                 ${PharaoDir}/getExt.cpp
                 ${PharaoDir}/hAccFuncCalc.cpp
                 ${PharaoDir}/hDonFuncCalc.cpp
                 ${PharaoDir}/hybridCalc.cpp
                 ${PharaoDir}/lipoFuncCalc.cpp
                 ${PharaoDir}/logOut.cpp
                 ${PharaoDir}/logPharmacophores.cpp
                 ${PharaoDir}/logScores.cpp
                 ${PharaoDir}/mainErr.cpp
                 ${PharaoDir}/mainWar.cpp
                 ${PharaoDir}/options.cpp
                 ${PharaoDir}/parseCommandLine.cpp
                 ${PharaoDir}/pharMerger.cpp
                 ${PharaoDir}/pharmacophore.cpp
                 ${PharaoDir}/printHeader.cpp
                 ${PharaoDir}/printInfo.cpp
                 ${PharaoDir}/printProgress.cpp
                 ${PharaoDir}/printUsage.cpp
                 ${PharaoDir}/result.cpp
                 ${PharaoDir}/siMath.cpp
                 ${PharaoDir}/solutionInfo.cpp
                 ${PharaoDir}/stringTokenizer.cpp
                 ${PharaoDir}/utilities.cpp
)

# Include files from pharao project
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${OPENBABEL2_INCLUDE_DIR})

# Executables
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} ${OPENBABEL2_LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
