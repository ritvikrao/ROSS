# Charades' dragonfly model (CODES-derived), matching it draw for draw;
# --gen-time is its weak-scaling form. See dragonfly.main.c.
add_executable(dragonfly dragonfly.main.c)
target_link_libraries(dragonfly PRIVATE ROSS m)
ROSS_TEST_SCHEDULERS(dragonfly)
install(TARGETS dragonfly RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
