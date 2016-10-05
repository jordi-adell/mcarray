# - Try to find Verbio
# Once done this will define
#  WIPP_FOUND - System has Verbio
#  WIPP_INCLUDE_DIRS - The Verbio include directories
#  WIPP_LIBRARIES - The libraries needed to use Verbio


find_path(WIPP_INCLUDE_DIR wipp/wipp.h)
find_library(WIPP_LIBRARY wipp)

list(APPEND WIPP_REQUIRED_VARS WIPP_INCLUDE_DIR)
list(APPEND WIPP_REQUIRED_VARS WIPP_LIBRARY)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(WIPP
				  REQUIRED_VARS ${WIPP_REQUIRED_VARS}
  )

mark_as_advanced(${WIPP_REQUIRED_VARS})

