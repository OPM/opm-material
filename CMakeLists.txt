# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# set up project
project("opm-material" CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.12)

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

# include the dune macros
include(DuneMacros)

# # include the OPM cmake macros
# include(OpmMacros)

# # start a dune project with information from dune.module
# dune_project()

# add_subdirectory("dune")
# add_subdirectory("doc")
# add_subdirectory("lib")
# add_subdirectory("src")
# add_subdirectory("cmake/modules")

# # finalize the dune project, e.g. generating config.h etc.
# finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

# opm_add_test(test_eclblackoilpvt CONDITION OPM_PARSER_FOUND)
# opm_add_test(test_eclmateriallawmanager CONDITION OPM_PARSER_FOUND)
# opm_add_test(test_fluidmatrixinteractions)
# opm_add_test(test_pengrobinson)
# opm_add_test(test_densead)
# opm_add_test(test_ncpflash)
# opm_add_test(test_spline)
# opm_add_test(test_tabulation)
# opm_add_test(test_2dtables)
# opm_add_test(test_components)
# opm_add_test(test_fluidsystems)
# opm_add_test(test_immiscibleflash)
