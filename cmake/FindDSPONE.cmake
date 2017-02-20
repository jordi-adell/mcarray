# - Try to find Verbio
# Once done this will define
#  DSPONE_FOUND - System has Verbio
#  DSPONE_INCLUDE_DIRS - The Verbio include directories
#  DSPONE_LIBRARIES - The libraries needed to use Verbio

set(CMAKE_FIND_ROOT_PATH ${CMAKE_INSTALL_PREFIX})

find_path(DSPONE_INCLUDE_DIRS dspone/dsp.h)
find_library(DSPONE_LIBRARIES dspone)

list(APPEND DSPONE_REQUIRED_VARS DSPONE_INCLUDE_DIRS)
list(APPEND DSPONE_REQUIRED_VARS DSPONE_LIBRARIES)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DSPONE
				  REQUIRED_VARS ${DSPONE_REQUIRED_VARS}
  )

mark_as_advanced(${DSPONE_REQUIRED_VARS})

find_path(GUI_H_DIR dspone/plot/dspgui.h)
if (GUI_H_DIR)
    set(DSPONE_GUI "DSPONE_GUI-FOUND")
else(GUI_H_DIR)
    set(DSPONE_GUI "DSPONE_GUI-NOTFOUND")
endif(GUI_H_DIR)
