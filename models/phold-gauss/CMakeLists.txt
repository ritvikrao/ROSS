# Gaussian-destination PHOLD (Bachan et al., PADS'24 / Devastator), matching
# Charades' models/phold/simple semantics. See phold-gauss.main.c.
add_executable(phold-gauss phold-gauss.main.c)
target_link_libraries(phold-gauss PRIVATE ROSS m)
ROSS_TEST_SCHEDULERS(phold-gauss)
install(TARGETS phold-gauss RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
