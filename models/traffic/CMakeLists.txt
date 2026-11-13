# Charades' traffic model (grid of intersections, cars), matching it draw for
# draw; --trip-radius is its weak-scaling form. See traffic.main.c.
add_executable(traffic traffic.main.c)
target_link_libraries(traffic PRIVATE ROSS m)
ROSS_TEST_SCHEDULERS(traffic)
install(TARGETS traffic RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
