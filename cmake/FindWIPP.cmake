# - Try to find Verbio
# Once done this will define
#  WIPP_FOUND - System has Verbio
#  WIPP_INCLUDE_DIRS - The Verbio include directories
#  WIPP_LIBRARIES - The libraries needed to use Verbio

set(CMAKE_FIND_ROOT_PATH ${CMAKE_INSTALL_PREFIX})

find_path(WIPP_INCLUDE_DIRS wipp/wipp.h)
find_library(WIPP_LIBRARIES wipp)

list(APPEND WIPP_REQUIRED_VARS WIPP_INCLUDE_DIRS)
list(APPEND WIPP_REQUIRED_VARS WIPP_LIBRARIES)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(WIPP
				  REQUIRED_VARS ${WIPP_REQUIRED_VARS}
  )

mark_as_advanced(${WIPP_REQUIRED_VARS})

