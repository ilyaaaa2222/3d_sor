module subgrid_solver
   implicit none
   private

   interface
      subroutine i_subgrid_solver_method(subgrid0, subgrid_size0)
         implicit none
         integer, intent(in) :: subgrid_size0
         real*8, intent(inout) :: subgrid0(subgrid_size0*subgrid_size0)
      end subroutine i_subgrid_solver_method
   end interface

   public :: i_subgrid_solver_method
   public :: set_subgrid_solver_settings, set_default_subgrid_solver_settings
   public :: get_subgrid_solver_settings, get_subgrid_solver_results
   public :: run_subgrid_solver

   public :: original_sor
   public :: tiling_sor
   public :: tiling_sor_1
   public :: tiling_sor_2
   public :: tiling_sor_4
   public :: tiling_sor_8
   public :: tiling_sor_16
   public :: subtiling_sor
   public :: subtiling_sor_test_version
   public :: subtiling_sor_1
   public :: subtiling_sor_2
   public :: subtiling_sor_4
   public :: subtiling_sor_8
   public :: subtiling_sor_16

   real*8 :: eps, omega
   integer :: max_iter
   integer :: tile_size, subtile_level

   real*8 :: time
   integer :: iter

contains

   subroutine set_default_subgrid_solver_settings()
      implicit none

      eps = 1.0d-8
      max_iter = 100000
      omega = -1
      tile_size = -1
      subtile_level = -1

   end subroutine set_default_subgrid_solver_settings

   subroutine set_subgrid_solver_settings(new_eps, new_max_iter, new_omega, new_tile_size, new_subtile_level)
     implicit none
     real*8, intent(in), optional :: new_eps, new_omega
     integer, intent(in), optional :: new_max_iter, new_tile_size, new_subtile_level
 
     if (present(new_eps)) eps = new_eps
     if (present(new_max_iter)) max_iter = new_max_iter
     if (present(new_omega)) omega = new_omega
     if (present(new_tile_size)) tile_size = new_tile_size
     if (present(new_subtile_level)) subtile_level = new_subtile_level
 
   end subroutine set_subgrid_solver_settings

   subroutine get_subgrid_solver_settings(new_eps, new_max_iter, new_omega, new_tile_size, new_subtile_level)
      implicit none
      real*8, intent(out), optional :: new_eps, new_omega
      integer, intent(out), optional :: new_max_iter, new_tile_size, new_subtile_level
  
      if (present(new_eps)) new_eps = eps
      if (present(new_max_iter)) new_max_iter = max_iter
      if (present(new_omega)) new_omega = omega
      if (present(new_tile_size)) new_tile_size = tile_size
      if (present(new_subtile_level)) new_subtile_level = subtile_level
  
    end subroutine get_subgrid_solver_settings
 
   subroutine run_subgrid_solver(new_subgrid, new_subgrid_size, new_subgrid_solver_method)
      implicit none
      integer, intent(in) :: new_subgrid_size
      real*8, intent(inout), target :: new_subgrid(:)
      procedure(i_subgrid_solver_method) :: new_subgrid_solver_method

      real*8 :: start_time, end_time

      call cpu_time(start_time)
      call new_subgrid_solver_method(new_subgrid, new_subgrid_size)
      call cpu_time(end_time)

      time = end_time - start_time

   end subroutine run_subgrid_solver

   subroutine get_subgrid_solver_results(res_time, res_iter)
      implicit none
      real*8, intent(out) :: res_time
      integer, intent(out) :: res_iter
      res_time = time
      res_iter = iter
   end subroutine get_subgrid_solver_results

   subroutine original_sor(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer :: i, l0, l1
      real*8 :: error, f0, f1, u_old

      error = eps + 1.0d0
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l1 = 3, u_size
            do l0 = 3, u_size

               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

               error = error + abs(u(i) - u_old)
               i = i + 1

            end do
            i = i + 2
         end do
         iter = iter + 1
      end do

   end subroutine original_sor

   subroutine tiling_sor(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer :: i, l0, l1, l2, l3, tile_count
      real*8 :: error, f0, f1, u_old

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size
               end do
               i = i - tile_size*(u_size - 1)
            end do
            i = i + u_size*(tile_size - 1) + 2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_sor

   subroutine tiling_sor_1(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 1
      integer :: i, l0, l1, l2, l3, tile_count
      real*8 :: error, f0, f1, u_old

      if (u_size <= tile_size_param) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size_param
               end do
               i = i - tile_size_param*(u_size - 1)
            end do
            i = i + u_size*(tile_size_param - 1) + 2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_sor_1

   subroutine tiling_sor_2(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 2
      integer :: i, l0, l1, l2, l3, tile_count
      real*8 :: error, f0, f1, u_old

      if (u_size <= tile_size_param) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size_param
               end do
               i = i - tile_size_param*(u_size - 1)
            end do
            i = i + u_size*(tile_size_param - 1) + 2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_sor_2

   subroutine tiling_sor_4(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 4
      integer :: i, l0, l1, l2, l3, tile_count
      real*8 :: error, f0, f1, u_old

      if (u_size <= tile_size_param) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size_param
               end do
               i = i - tile_size_param*(u_size - 1)
            end do
            i = i + u_size*(tile_size_param - 1) + 2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_sor_4

   subroutine tiling_sor_8(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 8
      integer :: i, l0, l1, l2, l3, tile_count
      real*8 :: error, f0, f1, u_old

      if (u_size <= tile_size_param) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size_param
               end do
               i = i - tile_size_param*(u_size - 1)
            end do
            i = i + u_size*(tile_size_param - 1) + 2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_sor_8

   subroutine tiling_sor_16(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 16
      integer :: i, l0, l1, l2, l3, tile_count
      real*8 :: error, f0, f1, u_old

      if (u_size <= tile_size_param) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size + 2
         error = 0.0d0

         do l3 = 1, tile_count
            do l2 = 1, tile_count
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - tile_size_param
               end do
               i = i - tile_size_param*(u_size - 1)
            end do
            i = i + u_size*(tile_size_param - 1) + 2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_sor_16

   subroutine subtiling_sor(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer :: tile_count, i, l0, l1, l2, l3, l4
      real*8 :: f0, f1, u_old, error

      tile_count = (u_size - 2)/tile_size

      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size + 2

         do l2 = 0, subtile_level
            do l1 = 1, tile_size - l2
               do l0 = 1, tile_size - l2

                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*(l1 - 1)
         end do
         i = i + tile_size

         do l3 = 2, tile_count - 1
            do l2 = 0, subtile_level
               do l1 = 1, tile_size - l2
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*(l1 - 1) - 1
            end do
            i = i + tile_size + subtile_level + 1
         end do

         do l2 = 0, subtile_level
            do l1 = 1, tile_size - l2
               do l0 = 1, tile_size + l2

                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*(l1 - 1) - 1
         end do
         i = i + tile_size + subtile_level + u_size*(tile_size - 1) + 3

         do l4 = 2, tile_count - 1
            do l2 = 0, subtile_level
               do l1 = 1, tile_size
                  do l0 = 1, tile_size - l2

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1
            end do
            i = i + u_size*(subtile_level + 1) + tile_size

            do l3 = 2, tile_count - 1
               do l2 = 0, subtile_level
                  do l1 = 1, tile_size
                     do l0 = 1, tile_size

                        u_old = u(i)
                        u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                        if (l2.eq.subtile_level) then
                           error = error + abs(u(i) - u_old)
                        end if
                        i = i + 1

                     end do
                     i = i + u_size - l0 + 1
                  end do
                  i = i - u_size*l1 - 1
               end do
               i = i + u_size*(subtile_level + 1) + tile_size + subtile_level + 1
            end do

            do l2 = 0, subtile_level
               do l1 = 1, tile_size
                  do l0 = 1, tile_size + l2

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(subtile_level + 1) + tile_size + subtile_level + u_size*(tile_size - 1) + 3
         end do

         do l2 = 0, subtile_level
            do l1 = 1, tile_size + l2
               do l0 = 1, tile_size - l2

                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1
         end do
         i = i + u_size*(subtile_level + 1) + tile_size

         do l3 = 2, tile_count - 1
            do l2 = 0, subtile_level
               do l1 = 1, tile_size + l2
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(subtile_level + 1) + tile_size + subtile_level + 1
         end do

         do l2 = 0, subtile_level
            do l1 = 1, tile_size + l2
               do l0 = 1, tile_size + l2

                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                  if (l2.eq.subtile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1 - 1
         end do

         if (error < eps) then
            return
         end if

      end do
   end subroutine subtiling_sor

   subroutine subtiling_sor_test_version(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, dimension(:,:), allocatable :: subtile_indices
      integer :: tile_count, subtile_count, subtile_idx, i, l0, l1, l2, l3, l4, s0, s1, s2, s3
      real*8 :: f0, f1, u_old, error

      tile_count = (u_size - 2)/tile_size
      subtile_count = tile_count*tile_count*(subtile_level+1)
      allocate(subtile_indices(subtile_count, 4))
      subtile_idx = 1

      do l4 = 1, tile_count

         if (l4.eq.1) then
            s1 = -1
         else if (l4.eq.tile_count) then
            s1 = 1
         else
            s1 = 0
         end if

         if (l4.eq.1) then
            s3 = 0
         else
            s3 = 1
         end if

         do l3 = 1, tile_count

            if (l3.eq.1) then
               s0 = -1
            else if (l3.eq.tile_count) then
               s0 = 1
            else
               s0 = 0
            end if

            if (l3.eq.1) then
               s2 = 0
            else
               s2 = 1
            end if

            do l2 = 0, subtile_level
               subtile_indices(subtile_idx, 1) = (l4-1)*tile_size*u_size + (l3-1)*tile_size + u_size + 2 - l2*s3*u_size - l2*s2
               subtile_indices(subtile_idx, 2) = tile_size + l2*s1
               subtile_indices(subtile_idx, 3) = tile_size + l2*s0
               if (l2.eq.subtile_level) then
                  subtile_indices(subtile_idx, 4) = 1
               else
                  subtile_indices(subtile_idx, 4) = 0
               end if
               subtile_idx = subtile_idx + 1

            end do
         end do
      end do

      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size + 2

         do subtile_idx = 1, subtile_count
            i = subtile_indices(subtile_idx, 1)
            l2 = subtile_indices(subtile_idx, 2)
            l3 = subtile_indices(subtile_idx, 3)
            l4 = subtile_indices(subtile_idx, 4)
            if (l4.eq.1) then
               do l1 = 1, l2
                  do l0 = 1, l3

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     error = error + abs(u(i) - u_old)
                     i = i + 1

                  end do
                  i = i + u_size - l3
               end do
            else
               do l1 = 1, l2
                  do l0 = 1, l3

                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))

                     i = i + 1

                  end do
                  i = i + u_size - l3
               end do
            end if
         end do

         if (error < eps) then
            exit
         end if

      end do

      deallocate(subtile_indices)

   end subroutine subtiling_sor_test_version

   subroutine subtiling_sor_1(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 1
      integer, parameter :: subtile_level_param = 1
      integer :: i, l0, l1, l2, l3, l4, tile_count
      real*8 :: error, f0, f1, u_old

      iter = 0
      if (u_size <= tile_size_param*2) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0
      
      do while (error > eps .and. iter < max_iter)
         i = u_size + 2
         error = 0.0d0
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i + tile_size_param - u_size*(tile_size_param - 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - tile_size_param*u_size - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(tile_size_param - 1)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i + subtile_level_param + tile_size_param + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 1
            end do
            i = i + subtile_level_param*u_size - tile_size_param*u_size + tile_size_param
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i + subtile_level_param*u_size + subtile_level_param - tile_size_param*u_size + tile_size_param
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 1
            end do
            i = i + subtile_level_param*u_size + subtile_level_param + tile_size_param - u_size + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1)
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i + subtile_level_param*u_size - tile_size_param*u_size + tile_size_param - u_size
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param*u_size + subtile_level_param - tile_size_param*u_size + tile_size_param - u_size
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1) - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param + 2) - 1
         iter = iter + subtile_level_param + 1
      end do

   end subroutine subtiling_sor_1

   subroutine subtiling_sor_2(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 2
      integer, parameter :: subtile_level_param = 2
      integer :: i, l0, l1, l2, l3, l4, tile_count
      real*8 :: error, f0, f1, u_old

      iter = 0
      if (u_size <= tile_size_param*2) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      do while (error > eps .and. iter < max_iter)
         i = u_size + 2
         error = 0.0d0
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i + tile_size_param - u_size*(tile_size_param - 2)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - tile_size_param*u_size - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(tile_size_param - 2)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i + subtile_level_param + tile_size_param + u_size + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 1
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 2
            end do
            i = i + subtile_level_param*u_size - tile_size_param*u_size + tile_size_param
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i + subtile_level_param*u_size + subtile_level_param - tile_size_param*u_size + tile_size_param
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 1
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 2
            end do
            i = i + subtile_level_param*u_size + subtile_level_param + tile_size_param - u_size + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1)
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param + 2)
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 3)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 2) - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 3)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1) - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param + 2) - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param + 3) - 1
         iter = iter + subtile_level_param + 1
      end do

   end subroutine subtiling_sor_2

   subroutine subtiling_sor_4(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 4
      integer, parameter :: subtile_level_param = 4
      integer :: i, l0, l1, l2, l3, l4, tile_count
      real*8 :: error, f0, f1, u_old

      iter = 0
      if (u_size <= tile_size_param*2) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      do while (error > eps .and. iter < max_iter)
         i = u_size + 2
         error = 0.0d0
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i - u_size*(tile_size_param - 2)
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 3
         end do
         i = i - u_size*(tile_size_param - 3)
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 4
         end do
         i = i + tile_size_param - u_size*(tile_size_param - 4)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - tile_size_param*u_size - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 2) - 1
            do l1 = 1, tile_size_param - 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 3) - 1
            do l1 = 1, tile_size_param - 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(tile_size_param - 4)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param - 2) - 1
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 3
         end do
         i = i - u_size*(tile_size_param - 3) - 1
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 4
         end do
         i = i + subtile_level_param + tile_size_param + 3*u_size + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 1
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 2
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 3
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 4
            end do
            i = i + subtile_level_param*u_size - tile_size_param*u_size + tile_size_param
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i + subtile_level_param*u_size + subtile_level_param - tile_size_param*u_size + tile_size_param
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 1
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 2
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 3
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 4
            end do
            i = i + subtile_level_param*u_size + subtile_level_param + tile_size_param - u_size + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1)
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param + 2)
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i - u_size*(tile_size_param + 3)
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 3
         end do
         i = i - u_size*(tile_size_param + 4)
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 4
         end do
         i = i + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 5)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 2) - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 3) - 1
            do l1 = 1, tile_size_param + 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 4) - 1
            do l1 = 1, tile_size_param + 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 5)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1) - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param + 2) - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param + 3) - 1
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 3
         end do
         i = i - u_size*(tile_size_param + 4) - 1
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 4
         end do
         i = i - u_size*(tile_size_param + 5) - 1
         iter = iter + subtile_level_param + 1
      end do

   end subroutine subtiling_sor_4

   subroutine subtiling_sor_8(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 8
      integer, parameter :: subtile_level_param = 8
      integer :: i, l0, l1, l2, l3, l4, tile_count
      real*8 :: error, f0, f1, u_old

      iter = 0
      if (u_size <= tile_size_param*2) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      do while (error > eps .and. iter < max_iter)
         i = u_size + 2
         error = 0.0d0
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i - u_size*(tile_size_param - 2)
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 3
         end do
         i = i - u_size*(tile_size_param - 3)
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 4
         end do
         i = i - u_size*(tile_size_param - 4)
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 5
         end do
         i = i - u_size*(tile_size_param - 5)
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 6
         end do
         i = i - u_size*(tile_size_param - 6)
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 7
         end do
         i = i - u_size*(tile_size_param - 7)
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 8
         end do
         i = i + tile_size_param - u_size*(tile_size_param - 8)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - tile_size_param*u_size - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 2) - 1
            do l1 = 1, tile_size_param - 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 3) - 1
            do l1 = 1, tile_size_param - 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 4) - 1
            do l1 = 1, tile_size_param - 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 5) - 1
            do l1 = 1, tile_size_param - 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 6) - 1
            do l1 = 1, tile_size_param - 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 7) - 1
            do l1 = 1, tile_size_param - 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(tile_size_param - 8)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param - 2) - 1
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 3
         end do
         i = i - u_size*(tile_size_param - 3) - 1
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 4
         end do
         i = i - u_size*(tile_size_param - 4) - 1
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 5
         end do
         i = i - u_size*(tile_size_param - 5) - 1
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 6
         end do
         i = i - u_size*(tile_size_param - 6) - 1
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 7
         end do
         i = i - u_size*(tile_size_param - 7) - 1
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 8
         end do
         i = i + subtile_level_param + tile_size_param + 7*u_size + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 1
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 2
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 3
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 4
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 5
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 6
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 7
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 8
            end do
            i = i + subtile_level_param*u_size - tile_size_param*u_size + tile_size_param
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i + subtile_level_param*u_size + subtile_level_param - tile_size_param*u_size + tile_size_param
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 1
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 2
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 3
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 4
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 5
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 6
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 7
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 8
            end do
            i = i + subtile_level_param*u_size + subtile_level_param + tile_size_param - u_size + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1)
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param + 2)
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i - u_size*(tile_size_param + 3)
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 3
         end do
         i = i - u_size*(tile_size_param + 4)
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 4
         end do
         i = i - u_size*(tile_size_param + 5)
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 5
         end do
         i = i - u_size*(tile_size_param + 6)
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 6
         end do
         i = i - u_size*(tile_size_param + 7)
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 7
         end do
         i = i - u_size*(tile_size_param + 8)
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 8
         end do
         i = i + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 9)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 2) - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 3) - 1
            do l1 = 1, tile_size_param + 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 4) - 1
            do l1 = 1, tile_size_param + 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 5) - 1
            do l1 = 1, tile_size_param + 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 6) - 1
            do l1 = 1, tile_size_param + 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 7) - 1
            do l1 = 1, tile_size_param + 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 8) - 1
            do l1 = 1, tile_size_param + 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 9)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1) - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param + 2) - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param + 3) - 1
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 3
         end do
         i = i - u_size*(tile_size_param + 4) - 1
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 4
         end do
         i = i - u_size*(tile_size_param + 5) - 1
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 5
         end do
         i = i - u_size*(tile_size_param + 6) - 1
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 6
         end do
         i = i - u_size*(tile_size_param + 7) - 1
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 7
         end do
         i = i - u_size*(tile_size_param + 8) - 1
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 8
         end do
         i = i - u_size*(tile_size_param + 9) - 1
         iter = iter + subtile_level_param + 1
      end do

   end subroutine subtiling_sor_8

   subroutine subtiling_sor_16(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size)

      integer, parameter :: tile_size_param = 16
      integer, parameter :: subtile_level_param = 16
      integer :: i, l0, l1, l2, l3, l4, tile_count
      real*8 :: error, f0, f1, u_old

      iter = 0
      if (u_size <= tile_size_param*2) return

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size_param
      f0 = 1.0d0 - omega
      f1 = omega * 0.25d0

      do while (error > eps .and. iter < max_iter)
         i = u_size + 2
         error = 0.0d0
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param - 1)
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i - u_size*(tile_size_param - 2)
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 3
         end do
         i = i - u_size*(tile_size_param - 3)
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 4
         end do
         i = i - u_size*(tile_size_param - 4)
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 5
         end do
         i = i - u_size*(tile_size_param - 5)
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 6
         end do
         i = i - u_size*(tile_size_param - 6)
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 7
         end do
         i = i - u_size*(tile_size_param - 7)
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 8
         end do
         i = i - u_size*(tile_size_param - 8)
         do l1 = 1, tile_size_param - 9
            do l0 = 1, tile_size_param - 9
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 9
         end do
         i = i - u_size*(tile_size_param - 9)
         do l1 = 1, tile_size_param - 10
            do l0 = 1, tile_size_param - 10
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 10
         end do
         i = i - u_size*(tile_size_param - 10)
         do l1 = 1, tile_size_param - 11
            do l0 = 1, tile_size_param - 11
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 11
         end do
         i = i - u_size*(tile_size_param - 11)
         do l1 = 1, tile_size_param - 12
            do l0 = 1, tile_size_param - 12
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 12
         end do
         i = i - u_size*(tile_size_param - 12)
         do l1 = 1, tile_size_param - 13
            do l0 = 1, tile_size_param - 13
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 13
         end do
         i = i - u_size*(tile_size_param - 13)
         do l1 = 1, tile_size_param - 14
            do l0 = 1, tile_size_param - 14
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 14
         end do
         i = i - u_size*(tile_size_param - 14)
         do l1 = 1, tile_size_param - 15
            do l0 = 1, tile_size_param - 15
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 15
         end do
         i = i - u_size*(tile_size_param - 15)
         do l1 = 1, tile_size_param - 16
            do l0 = 1, tile_size_param - 16
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 16
         end do
         i = i + tile_size_param - u_size*(tile_size_param - 16)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param - 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - tile_size_param*u_size - 1
            do l1 = 1, tile_size_param - 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 1) - 1
            do l1 = 1, tile_size_param - 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 2) - 1
            do l1 = 1, tile_size_param - 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 3) - 1
            do l1 = 1, tile_size_param - 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 4) - 1
            do l1 = 1, tile_size_param - 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 5) - 1
            do l1 = 1, tile_size_param - 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 6) - 1
            do l1 = 1, tile_size_param - 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 7) - 1
            do l1 = 1, tile_size_param - 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 8) - 1
            do l1 = 1, tile_size_param - 9
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 9) - 1
            do l1 = 1, tile_size_param - 10
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 10) - 1
            do l1 = 1, tile_size_param - 11
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 11) - 1
            do l1 = 1, tile_size_param - 12
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 12) - 1
            do l1 = 1, tile_size_param - 13
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 13) - 1
            do l1 = 1, tile_size_param - 14
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 14) - 1
            do l1 = 1, tile_size_param - 15
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param - 15) - 1
            do l1 = 1, tile_size_param - 16
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param - u_size*(tile_size_param - 16)
         end do
         do l1 = 1, tile_size_param - 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - tile_size_param*u_size - 1
         do l1 = 1, tile_size_param - 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param - 1) - 1
         do l1 = 1, tile_size_param - 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param - 2) - 1
         do l1 = 1, tile_size_param - 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 3
         end do
         i = i - u_size*(tile_size_param - 3) - 1
         do l1 = 1, tile_size_param - 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 4
         end do
         i = i - u_size*(tile_size_param - 4) - 1
         do l1 = 1, tile_size_param - 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 5
         end do
         i = i - u_size*(tile_size_param - 5) - 1
         do l1 = 1, tile_size_param - 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 6
         end do
         i = i - u_size*(tile_size_param - 6) - 1
         do l1 = 1, tile_size_param - 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 7
         end do
         i = i - u_size*(tile_size_param - 7) - 1
         do l1 = 1, tile_size_param - 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 8
         end do
         i = i - u_size*(tile_size_param - 8) - 1
         do l1 = 1, tile_size_param - 9
            do l0 = 1, tile_size_param + 9
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 9
         end do
         i = i - u_size*(tile_size_param - 9) - 1
         do l1 = 1, tile_size_param - 10
            do l0 = 1, tile_size_param + 10
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 10
         end do
         i = i - u_size*(tile_size_param - 10) - 1
         do l1 = 1, tile_size_param - 11
            do l0 = 1, tile_size_param + 11
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 11
         end do
         i = i - u_size*(tile_size_param - 11) - 1
         do l1 = 1, tile_size_param - 12
            do l0 = 1, tile_size_param + 12
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 12
         end do
         i = i - u_size*(tile_size_param - 12) - 1
         do l1 = 1, tile_size_param - 13
            do l0 = 1, tile_size_param + 13
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 13
         end do
         i = i - u_size*(tile_size_param - 13) - 1
         do l1 = 1, tile_size_param - 14
            do l0 = 1, tile_size_param + 14
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 14
         end do
         i = i - u_size*(tile_size_param - 14) - 1
         do l1 = 1, tile_size_param - 15
            do l0 = 1, tile_size_param + 15
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 15
         end do
         i = i - u_size*(tile_size_param - 15) - 1
         do l1 = 1, tile_size_param - 16
            do l0 = 1, tile_size_param + 16
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 16
         end do
         i = i + subtile_level_param + tile_size_param + 15*u_size + 2
         do l4 = 2, tile_count - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 1
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 2
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 3
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 4
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 5
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 6
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 7
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 8
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 9
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 9
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 10
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 10
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 11
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 11
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 12
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 12
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 13
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 13
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 14
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 14
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 15
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 15
            end do
            i = i - u_size*(tile_size_param + 1)
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param - 16
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size + 16
            end do
            i = i + subtile_level_param*u_size - tile_size_param*u_size + tile_size_param
            do l3 = 2, tile_count - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i - u_size*(tile_size_param + 1) - 1
               do l1 = 1, tile_size_param
                  do l0 = 1, tile_size_param
                     u_old = u(i)
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                     error = error + Abs(u_old - u(i))
                     i = i + 1
                  end do
                  i = i - tile_size_param + u_size
               end do
               i = i + subtile_level_param*u_size + subtile_level_param - tile_size_param*u_size + tile_size_param
            end do
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 0
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 1
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 1
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 2
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 2
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 3
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 3
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 4
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 4
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 5
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 5
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 6
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 6
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 7
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 7
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 8
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 8
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 9
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 9
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 10
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 10
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 11
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 11
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 12
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 12
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 13
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 13
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 14
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 14
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 15
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 15
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param
               do l0 = 1, tile_size_param + 16
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size - 16
            end do
            i = i + subtile_level_param*u_size + subtile_level_param + tile_size_param - u_size + 2
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1)
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 1
         end do
         i = i - u_size*(tile_size_param + 2)
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 2
         end do
         i = i - u_size*(tile_size_param + 3)
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 3
         end do
         i = i - u_size*(tile_size_param + 4)
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 4
         end do
         i = i - u_size*(tile_size_param + 5)
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 5
         end do
         i = i - u_size*(tile_size_param + 6)
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 6
         end do
         i = i - u_size*(tile_size_param + 7)
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 7
         end do
         i = i - u_size*(tile_size_param + 8)
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 8
         end do
         i = i - u_size*(tile_size_param + 9)
         do l1 = 1, tile_size_param + 9
            do l0 = 1, tile_size_param - 9
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 9
         end do
         i = i - u_size*(tile_size_param + 10)
         do l1 = 1, tile_size_param + 10
            do l0 = 1, tile_size_param - 10
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 10
         end do
         i = i - u_size*(tile_size_param + 11)
         do l1 = 1, tile_size_param + 11
            do l0 = 1, tile_size_param - 11
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 11
         end do
         i = i - u_size*(tile_size_param + 12)
         do l1 = 1, tile_size_param + 12
            do l0 = 1, tile_size_param - 12
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 12
         end do
         i = i - u_size*(tile_size_param + 13)
         do l1 = 1, tile_size_param + 13
            do l0 = 1, tile_size_param - 13
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 13
         end do
         i = i - u_size*(tile_size_param + 14)
         do l1 = 1, tile_size_param + 14
            do l0 = 1, tile_size_param - 14
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 14
         end do
         i = i - u_size*(tile_size_param + 15)
         do l1 = 1, tile_size_param + 15
            do l0 = 1, tile_size_param - 15
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 15
         end do
         i = i - u_size*(tile_size_param + 16)
         do l1 = 1, tile_size_param + 16
            do l0 = 1, tile_size_param - 16
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size + 16
         end do
         i = i + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 17)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size_param + 0
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 1) - 1
            do l1 = 1, tile_size_param + 1
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 2) - 1
            do l1 = 1, tile_size_param + 2
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 3) - 1
            do l1 = 1, tile_size_param + 3
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 4) - 1
            do l1 = 1, tile_size_param + 4
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 5) - 1
            do l1 = 1, tile_size_param + 5
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 6) - 1
            do l1 = 1, tile_size_param + 6
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 7) - 1
            do l1 = 1, tile_size_param + 7
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 8) - 1
            do l1 = 1, tile_size_param + 8
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 9) - 1
            do l1 = 1, tile_size_param + 9
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 10) - 1
            do l1 = 1, tile_size_param + 10
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 11) - 1
            do l1 = 1, tile_size_param + 11
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 12) - 1
            do l1 = 1, tile_size_param + 12
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 13) - 1
            do l1 = 1, tile_size_param + 13
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 14) - 1
            do l1 = 1, tile_size_param + 14
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 15) - 1
            do l1 = 1, tile_size_param + 15
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i - u_size*(tile_size_param + 16) - 1
            do l1 = 1, tile_size_param + 16
               do l0 = 1, tile_size_param
                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - tile_size_param + u_size
            end do
            i = i + subtile_level_param + tile_size_param + u_size*(subtile_level_param + 1) - u_size*(tile_size_param + 17)
         end do
         do l1 = 1, tile_size_param + 0
            do l0 = 1, tile_size_param + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size
         end do
         i = i - u_size*(tile_size_param + 1) - 1
         do l1 = 1, tile_size_param + 1
            do l0 = 1, tile_size_param + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 1
         end do
         i = i - u_size*(tile_size_param + 2) - 1
         do l1 = 1, tile_size_param + 2
            do l0 = 1, tile_size_param + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 2
         end do
         i = i - u_size*(tile_size_param + 3) - 1
         do l1 = 1, tile_size_param + 3
            do l0 = 1, tile_size_param + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 3
         end do
         i = i - u_size*(tile_size_param + 4) - 1
         do l1 = 1, tile_size_param + 4
            do l0 = 1, tile_size_param + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 4
         end do
         i = i - u_size*(tile_size_param + 5) - 1
         do l1 = 1, tile_size_param + 5
            do l0 = 1, tile_size_param + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 5
         end do
         i = i - u_size*(tile_size_param + 6) - 1
         do l1 = 1, tile_size_param + 6
            do l0 = 1, tile_size_param + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 6
         end do
         i = i - u_size*(tile_size_param + 7) - 1
         do l1 = 1, tile_size_param + 7
            do l0 = 1, tile_size_param + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 7
         end do
         i = i - u_size*(tile_size_param + 8) - 1
         do l1 = 1, tile_size_param + 8
            do l0 = 1, tile_size_param + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 8
         end do
         i = i - u_size*(tile_size_param + 9) - 1
         do l1 = 1, tile_size_param + 9
            do l0 = 1, tile_size_param + 9
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 9
         end do
         i = i - u_size*(tile_size_param + 10) - 1
         do l1 = 1, tile_size_param + 10
            do l0 = 1, tile_size_param + 10
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 10
         end do
         i = i - u_size*(tile_size_param + 11) - 1
         do l1 = 1, tile_size_param + 11
            do l0 = 1, tile_size_param + 11
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 11
         end do
         i = i - u_size*(tile_size_param + 12) - 1
         do l1 = 1, tile_size_param + 12
            do l0 = 1, tile_size_param + 12
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 12
         end do
         i = i - u_size*(tile_size_param + 13) - 1
         do l1 = 1, tile_size_param + 13
            do l0 = 1, tile_size_param + 13
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 13
         end do
         i = i - u_size*(tile_size_param + 14) - 1
         do l1 = 1, tile_size_param + 14
            do l0 = 1, tile_size_param + 14
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 14
         end do
         i = i - u_size*(tile_size_param + 15) - 1
         do l1 = 1, tile_size_param + 15
            do l0 = 1, tile_size_param + 15
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 15
         end do
         i = i - u_size*(tile_size_param + 16) - 1
         do l1 = 1, tile_size_param + 16
            do l0 = 1, tile_size_param + 16
               u_old = u(i)
               u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size))
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - tile_size_param + u_size - 16
         end do
         i = i - u_size*(tile_size_param + 17) - 1
         iter = iter + subtile_level_param + 1
      end do

   end subroutine subtiling_sor_16

end module subgrid_solver
