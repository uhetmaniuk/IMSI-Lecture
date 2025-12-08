//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <Kokkos_Core.hpp>

// Functors for CUDA nested parallelism (must be at namespace scope)
struct CudaTeamFunctor3 {
  Kokkos::View<int, Kokkos::CudaSpace> count;
  typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;

  KOKKOS_FUNCTION
  void operator()(const team_member& thread) const {
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(thread, 3),
        [=](const int& i) { Kokkos::atomic_fetch_add(&count(), 1); });
  }
};

struct CudaTeamFunctor5 {
  Kokkos::View<int, Kokkos::CudaSpace> count;
  typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;

  KOKKOS_FUNCTION
  void operator()(const team_member& thread) const {
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(thread, 5),
        [=](const int& i) { Kokkos::atomic_fetch_add(&count(), 1); });
  }
};

int
main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    {
      /// OpenMP - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, Kokkos::AUTO);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int, Kokkos::HostSpace>                                       count("Count");
      Kokkos::parallel_for(
          "OpenMP1", policy, KOKKOS_LAMBDA(const team_member& thread) { Kokkos::atomic_fetch_add(&count(), 1); });
      printf(
          " %d -- 1 parallel_for (27, AUTO) -- Kokkos Sizes%d %d  \n",
          count(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// OpenMP - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, 3);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int, Kokkos::HostSpace>                                       count("Count");
      Kokkos::parallel_for(
          "OpenMP2", policy, KOKKOS_LAMBDA(const team_member& thread) { Kokkos::atomic_fetch_add(&count(), 1); });
      printf(
          " %d -- 1 parallel_for (27, 3) -- Kokkos Sizes %d %d  \n", count(), policy.league_size(), policy.team_size());
    }

    {
      /// OpenMP - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, Kokkos::AUTO);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int, Kokkos::HostSpace>                                       count("Count");
      Kokkos::parallel_for(
          "OpenMP3", policy, KOKKOS_LAMBDA(const team_member& thread) {
            for (int k = 0; k < 3; ++k) { Kokkos::atomic_fetch_add(&count(), 1); }
          });
      printf(
          "OpenMP parallel_for - for (27, Auto) & 3 %d -- Kokkos Sizes %d %d  \n",
          count(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// OpenMP - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, 3);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int, Kokkos::HostSpace>                                       count("Count");
      Kokkos::parallel_for(
          "OpenMP4", policy, KOKKOS_LAMBDA(const team_member& thread) {
            for (int k = 0; k < 5; ++k) { Kokkos::atomic_fetch_add(&count(), 1); }
          });
      printf(
          "OpenMP parallel_for - for (27, 3) & 5 %d -- Kokkos Sizes %d %d  \n",
          count(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// OpenMP - TeamPolicy with functor to avoid nested lambda
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, Kokkos::AUTO);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int>                                       count("Count");

      // Functor for outer parallel_for
      struct TeamFunctor {
        Kokkos::View<int> count;

        KOKKOS_FUNCTION
        void operator()(const team_member& thread) const {
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(thread, 3),
              [=](const int& i) { Kokkos::atomic_fetch_add(&count(), 1); });
        }
      };

      Kokkos::parallel_for("OpenMP5", policy, TeamFunctor{count});
      printf(
          "OpenMP 2 parallel_for (27, Auto) & 3 %d -- Kokkos Sizes %d %d  \n",
          count(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// OpenMP - TeamPolicy with functor to avoid nested lambda
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, 3);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int>                                       count("Count");

      // Functor for outer parallel_for
      struct TeamFunctor {
        Kokkos::View<int> count;

        KOKKOS_FUNCTION
        void operator()(const team_member& thread) const {
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(thread, 5),
              [=](const int& i) { Kokkos::atomic_fetch_add(&count(), 1); });
        }
      };

      Kokkos::parallel_for("OpenMP6", policy, TeamFunctor{count});
      printf(
          "OpenMP 2 parallel_for (27, 3) & 5 %d -- Kokkos Sizes %d %d  \n",
          count(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// OpenMP - TeamPolicy with functor to avoid nested lambda
      Kokkos::TeamPolicy<Kokkos::OpenMP>                      policy(27, 7);
      typedef Kokkos::TeamPolicy<Kokkos::OpenMP>::member_type team_member;
      Kokkos::View<int>                                       count("Count");

      // Functor for outer parallel_for
      struct TeamFunctor {
        Kokkos::View<int> count;

        KOKKOS_FUNCTION
        void operator()(const team_member& thread) const {
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(thread, 5),
              [=](const int& i) { Kokkos::atomic_fetch_add(&count(), 1); });
        }
      };

      Kokkos::parallel_for("OpenMP7", policy, TeamFunctor{count});
      printf(
          "OpenMP 2 parallel_for (27, 7) & 5 %d -- Kokkos Sizes %d %d  \n",
          count(),
          policy.league_size(),
          policy.team_size());
    }
#ifdef KOKKOS_ENABLE_CUDA
    {
      /// CUDA - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, Kokkos::AUTO);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);
      Kokkos::parallel_for(
          "CUDA1", policy, KOKKOS_LAMBDA(const team_member& thread) { Kokkos::atomic_fetch_add(&count(), 1); });
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          " %d -- 1 parallel_for (27, AUTO) -- Kokkos Sizes%d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// CUDA - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, 3);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);
      Kokkos::parallel_for(
          "CUDA2", policy, KOKKOS_LAMBDA(const team_member& thread) { Kokkos::atomic_fetch_add(&count(), 1); });
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          " %d -- 1 parallel_for (27, 3) -- Kokkos Sizes %d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// CUDA - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, Kokkos::AUTO);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);
      Kokkos::parallel_for(
          "CUDA1", policy, KOKKOS_LAMBDA(const team_member& thread) {
            for (int k = 0; k < 3; ++k) { Kokkos::atomic_fetch_add(&count(), 1); }
          });
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          "CUDA parallel_for - for (27, Auto) & 3 %d -- Kokkos Sizes %d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// CUDA - TeamPolicy
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, 3);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);
      Kokkos::parallel_for(
          "CUDA1", policy, KOKKOS_LAMBDA(const team_member& thread) {
            for (int k = 0; k < 5; ++k) { Kokkos::atomic_fetch_add(&count(), 1); }
          });
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          "CUDA parallel_for - for (27, 3) & 5 %d -- Kokkos Sizes %d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// CUDA - TeamPolicy with functor to avoid nested lambda
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, Kokkos::AUTO);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);

      Kokkos::parallel_for("CUDA5", policy, CudaTeamFunctor3{count});
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          "CUDA 2 parallel_for (27, Auto) & 3 %d -- Kokkos Sizes %d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// CUDA - TeamPolicy with functor to avoid nested lambda
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, 3);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);

      Kokkos::parallel_for("CUDA6", policy, CudaTeamFunctor5{count});
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          "CUDA 2 parallel_for (27, 3) & 5 %d -- Kokkos Sizes %d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }

    {
      /// CUDA - TeamPolicy with functor to avoid nested lambda
      Kokkos::TeamPolicy<Kokkos::Cuda>                      policy(27, 7);
      typedef Kokkos::TeamPolicy<Kokkos::Cuda>::member_type team_member;
      Kokkos::View<int, Kokkos::CudaSpace>                  count("Count");
      auto                                                  count_host = Kokkos::create_mirror_view(count);

      Kokkos::parallel_for("CUDA7", policy, CudaTeamFunctor5{count});
      Kokkos::fence();
      Kokkos::deep_copy(count_host, count);
      printf(
          "CUDA 2 parallel_for (27, 7) & 5 %d -- Kokkos Sizes %d %d  \n",
          count_host(),
          policy.league_size(),
          policy.team_size());
    }
#endif
  }
  Kokkos::finalize();

  return 0;
}
