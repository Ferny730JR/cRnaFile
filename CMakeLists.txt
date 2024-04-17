cmake_minimum_required(VERSION 3.15.0)

project(rnaf C)

set(VERSION "0.1")

find_package(ZLIB REQUIRED)

option(RNAF_BUILD_EXAMPLES "Enable rnaf examples" OFF)

set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Installation directory for headers")
set(INSTALL_MAN_DIR "${CMAKE_INSTALL_PREFIX}/share/man" CACHE PATH "Installation directory for manual pages")

include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckCSourceCompiles)
enable_testing()

###################################################################################################
#  rnaf                                                                                           #
###################################################################################################

set(RNAF_PUBLIC_HEADERS
	include/rnaf.h)

set(RNAF_PRIVATE_HEADERS
	source/memory_utils.h
	source/string_utils.h)

set(RNAF_SOURCES
	source/rnaf.c
	source/memory_utils.c
	source/string_utils.c)

add_library(rnaf STATIC ${RNAF_SOURCES} ${RNAF_PUBLIC_HEADERS} ${RNAF_PRIVATE_HEADERS})
target_include_directories(rnaf PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} include)
target_link_libraries(rnaf ZLIB::ZLIB)

if(NOT SKIP_INSTALL_LIBRARIES AND NOT SKIP_INSTALL_ALL )
    install(TARGETS rnaf
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
        ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" )
endif()
if(NOT SKIP_INSTALL_HEADERS AND NOT SKIP_INSTALL_ALL )
    install(FILES ${RNAF_PUBLIC_HEADERS} DESTINATION "${INSTALL_INC_DIR}")
endif()

# Specify compiler options based on OS
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	target_compile_options(rnaf PRIVATE -Wall -O3)
	target_link_options(rnaf PRIVATE -Wall -O3 -flto)

elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin") # Apple
	target_compile_options(rnaf PRIVATE -Wall -Wextra -Wpedantic -O3)
	target_link_options(rnaf PRIVATE -Wall -Wextra -Wpedantic -O3 -flto)

else()
    message( SEND_ERROR "System ${CMAKE_SYSTEM_NAME} currently not supported.")
endif()
