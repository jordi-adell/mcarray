# - Try to find SNDFILE lirbary
# Once done this will define
#  SNDFILE_FOUND - System has SNDFILE
#  SNDFILE_INCLUDE_DIRS - The SNDFILE include directories
#  SNDFILE_LIBRARIES - The libraries needed to use SNDFILE


find_file(SNDFILE_INCLUDE_DIR sndfile.h)
find_library(SNDFILE_LIBRARIES sndfile)

list(APPEND SNDFILE_REQUIRED_VARS SNDFILE_INCLUDE_DIR)
list(APPEND SNDFILE_REQUIRED_VARS SNDFILE_LIBRARIES)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SNDFILE
				  REQUIRED_VARS ${SNDFILE_REQUIRED_VARS}
  )

mark_as_advanced(${SNDFILE_REQUIRED_VARS})
