add_subdirectory(kokkos)
if(USE_TACHO)
  add_subdirectory(tacho)
  message(STATUS "Tacho sparse solver library enabled")
else()
  message(STATUS "Tacho sparse solver library disabled")
endif()

add_subdirectory(kokkos-kernels)

### Disable for compilation times
#add_subdirectory(Vc)
######
