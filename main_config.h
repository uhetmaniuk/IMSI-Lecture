#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

#if defined(IMSI_USE_TACHO)
#include <Tacho.hpp>
#include <Tacho_CrsMatrixBase.hpp>
#include <Tacho_Driver.hpp>
#include <Tacho_Solver.hpp>
#endif

using accelerator_space = Kokkos::DefaultExecutionSpace;
using accelerator_type  = Kokkos::Device<accelerator_space, accelerator_space::memory_space>;

using host_execution_space = Kokkos::DefaultHostExecutionSpace;
using host_execution_type  = Kokkos::Device<host_execution_space, host_execution_space::memory_space>;

using serial_execution_space = Kokkos::Serial::execution_space;
using serial_execution_type  = Kokkos::Device<serial_execution_space, serial_execution_space::memory_space>;

