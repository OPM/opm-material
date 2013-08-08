# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# defines that must be present in config.h for our headers
set (opm-porsol_CONFIG_VAR
	HAVE_DYNAMIC_BOOST_TEST
	HAVE_VALGRIND
	)

# dependencies
set (opm-porsol_DEPS
	# compile with C99 support if available
	"C99"
	# compile with C++0x/11 support if available
	"CXX11Features"
	# various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
	# DUNE dependency
	"dune-common REQUIRED;
	dune-istl REQUIRED;
	dune-grid REQUIRED;
	opm-core REQUIRED;
	dune-cornerpoint REQUIRED"
	)
