# - Try to find Verbio
# Once done this will define
#  DSPONE_FOUND - System has Verbio
#  DSPONE_INCLUDE_DIRS - The Verbio include directories
#  DSPONE_LIBRARIES - The libraries needed to use Verbio


find_path(DSPONE_INCLUDE_DIR dspone/dsp.h)
find_library(DSPONE_LIBRARY dspone)

list(APPEND DSPONE_REQUIRED_VARS DSPONE_INCLUDE_DIR)
list(APPEND DSPONE_REQUIRED_VARS DSPONE_LIBRARY)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DSPONE
				  REQUIRED_VARS ${DSPONE_REQUIRED_VARS}
  )

mark_as_advanced(${DSPONE_REQUIRED_VARS})

