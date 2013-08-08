# - Turn on warnings when compiling

include (AddOptions)
include (UseCompVer)
is_compiler_gcc_compatible ()

if (CXX_COMPAT_GCC)
  # default warnings flags, if not set by user
  set_default_option (_warn_flag "-Wall" "(^|\ )-W")
  if (_warn_flag)
	message (STATUS "All warnings enabled: ${_warn_flag}")
	add_options (ALL_LANGUAGES ALL_BUILDS "${_warn_flag}")
  endif (_warn_flag)
endif ()
