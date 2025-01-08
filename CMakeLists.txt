cmake_minimum_required(VERSION 3.15)

project(
    Interpolation
    VERSION 0.0.1
    LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)    

include(FetchContent)
FetchContent_Declare(
  NumericConcepts
  GIT_REPOSITORY https://github.com/da380/NumericConcepts.git
  GIT_TAG main
)
FetchContent_MakeAvailable(NumericConcepts)

include(FetchContent)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
)
FetchContent_MakeAvailable(Eigen3)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(Interpolation INTERFACE NumericConcepts Eigen3::Eigen)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(BUILD_EXAMPLES   "whether or not examples should be built" ON)
    option(BUILD_TESTS      "whether or not tests should be built" ON)

    if(BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
    if(BUILD_TESTS)
        enable_testing()
        add_subdirectory(tests)
    endif()
endif()

include(GNUInstallDirs)

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)




