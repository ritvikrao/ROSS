# Charades' VLSI model (gate-level simulation of an ITC'99 netlist, copies
# joined in a ring). It compiles Charades' own netlist.C, so it needs the
# Charades tree: VLSI_NETLIST_DIR points at charades/models/vlsi. Skipped
# when that directory is missing. See vlsi.main.cpp.
set(VLSI_NETLIST_DIR "${CMAKE_SOURCE_DIR}/../../charades/models/vlsi"
    CACHE PATH "Charades models/vlsi (netlist.h, netlist.C)")
if(NOT EXISTS "${VLSI_NETLIST_DIR}/netlist.C")
  message(STATUS "vlsi model skipped: no netlist.C in VLSI_NETLIST_DIR (${VLSI_NETLIST_DIR})")
  return()
endif()
enable_language(CXX)
# netlist.C is C++ despite its extension (Charades' convention)
set_source_files_properties("${VLSI_NETLIST_DIR}/netlist.C" PROPERTIES LANGUAGE CXX)
add_executable(vlsi vlsi.main.cpp "${VLSI_NETLIST_DIR}/netlist.C")
target_include_directories(vlsi PRIVATE "${VLSI_NETLIST_DIR}")
set_target_properties(vlsi PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON LINKER_LANGUAGE CXX)
target_link_libraries(vlsi PRIVATE ROSS m)
install(TARGETS vlsi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
