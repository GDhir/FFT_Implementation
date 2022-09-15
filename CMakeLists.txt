cmake_minimum_required(VERSION 3.10)
Project(fft_impl)

set(${BOOST_ROOT} "/usr/include/boost")
set(${Boost_INCLUDE_DIR} "/usr/include/")
set(${CMAKE_BUILD_TYPE}, "Debug")

list (APPEND LIBS X11 GL dl)

set(myDesiredRootForHeaderFileSearch "/home/gaurav/Documents/DFT_Folder")
file(GLOB my_include_directories ${myDesiredRootForHeaderFileSearch})
include_directories(${my_include_directories})

# packages.
find_package(Boost REQUIRED COMPONENTS
    iostreams system filesystem
)

# Target.
add_library(gnuplot_iostream INTERFACE)
target_include_directories(gnuplot_iostream INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:>)
target_link_libraries(gnuplot_iostream INTERFACE
    Boost::iostreams
    Boost::system
    Boost::filesystem
)
add_executable(fft_impl fft_impl.cpp)
target_compile_features(fft_impl PRIVATE cxx_std_17)
target_compile_options(fft_impl PRIVATE -Wall -Wextra)
target_link_libraries(fft_impl PRIVATE gnuplot_iostream ${LIBS})