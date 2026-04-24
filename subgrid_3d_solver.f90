module subgrid_3d_solver
   implicit none
   private

   interface

      !Интерфейс метода решения задачи в трехмерной области.
      !subgrid0 - Массив значений в узлах кубической сетки.
      !subgrid_size0 - Размер кубической сетки вдоль одной оси.
      subroutine i_subgrid_3d_solver_method(subgrid0, subgrid_size0)
         implicit none
         integer, intent(in) :: subgrid_size0
         real*8, intent(inout), target :: subgrid0(subgrid_size0*subgrid_size0*subgrid_size0)
      end subroutine i_subgrid_3d_solver_method
   end interface

   !Методы с настройками решателя.
   public :: i_subgrid_3d_solver_method
   public :: set_subgrid_3d_solver_settings, set_default_subgrid_3d_solver_settings
   public :: get_subgrid_3d_solver_settings, get_subgrid_3d_solver_results
   public :: run_subgrid_3d_solver

   !Методы решения.
   public :: original_3d_sor
   public :: tiling_3d_sor
   public :: subtiling_3d_sor, subtiling_3d_sor_test

   !Параметры метода решения.
   real*8 :: eps, omega
   integer :: max_iter
   integer :: tile_size, subtile_level

   real*8 :: time
   integer :: iter

contains

   !Устанавливает настройки решателя по умолчанию.
   subroutine set_default_subgrid_3d_solver_settings()
      implicit none

      eps = 1.0d-8
      max_iter = 100000
      omega = -1
      tile_size = -1
      subtile_level = -1

   end subroutine set_default_subgrid_3d_solver_settings

   !Устанавливает настройки решателя.
   subroutine set_subgrid_3d_solver_settings(new_eps, new_max_iter, new_omega, new_tile_size, new_subtile_level)
     implicit none
     real*8, intent(in), optional :: new_eps, new_omega
     integer, intent(in), optional :: new_max_iter, new_tile_size, new_subtile_level
 
     if (present(new_eps)) eps = new_eps
     if (present(new_max_iter)) max_iter = new_max_iter
     if (present(new_omega)) omega = new_omega
     if (present(new_tile_size)) tile_size = new_tile_size
     if (present(new_subtile_level)) subtile_level = new_subtile_level
 
   end subroutine set_subgrid_3d_solver_settings

   !Возвращает параметры решателя.
   subroutine get_subgrid_3d_solver_settings(new_eps, new_max_iter, new_omega, new_tile_size, new_subtile_level)
      implicit none
      real*8, intent(out), optional :: new_eps, new_omega
      integer, intent(out), optional :: new_max_iter, new_tile_size, new_subtile_level
  
      if (present(new_eps)) new_eps = eps
      if (present(new_max_iter)) new_max_iter = max_iter
      if (present(new_omega)) new_omega = omega
      if (present(new_tile_size)) new_tile_size = tile_size
      if (present(new_subtile_level)) new_subtile_level = subtile_level
  
   end subroutine get_subgrid_3d_solver_settings
 
   !Запускает решатель.
   !new_subgrid - значения в узлах сетки.
   !new_subgrid_size - число узлов сетки вдоль одной оси.
   subroutine run_subgrid_3d_solver(new_subgrid, new_subgrid_size, new_subgrid_solver_method)
      implicit none
      integer, intent(in) :: new_subgrid_size
      real*8, intent(inout), target :: new_subgrid(:)
      procedure(i_subgrid_3d_solver_method) :: new_subgrid_solver_method

      real*8 :: start_time, end_time

      call cpu_time(start_time)
      call new_subgrid_solver_method(new_subgrid, new_subgrid_size)
      call cpu_time(end_time)

      time = end_time - start_time

   end subroutine run_subgrid_3d_solver

   !Возвращает результаты решателя.
   subroutine get_subgrid_3d_solver_results(res_time, res_iter)
      implicit none
      real*8, intent(out) :: res_time
      integer, intent(out) :: res_iter
      res_time = time
      res_iter = iter
   end subroutine get_subgrid_3d_solver_results

   !Метод последовательной верхней релаксации
   !с классической схемой обхода узлов сетки.
   subroutine original_3d_sor(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout), target :: u(u_size*u_size*u_size)

      integer :: i, l0, l1, l2, u_size2
      real*8 :: error, f0, f1, u_old

      u_size2 = u_size*u_size

      error = eps + 1.0d0
      f0 = 1.0d0 - omega
      f1 = omega / 6.0d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size*(u_size + 1) + 2
         error = 0.0d0

         do l2 = 3, u_size
            do l1 = 3, u_size
               do l0 = 3, u_size

                  u_old = u(i)
                  u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                     + u(i+u_size2) + u(i-u_size2))

                  error = error + abs(u(i) - u_old)
                  i = i + 1

               end do
               i = i + 2
            end do
            i = i + u_size*2
         end do
         iter = iter + 1
      end do

   end subroutine original_3d_sor

   !SOR с тайлингом.
   subroutine tiling_3d_sor(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout), target :: u(u_size*u_size*u_size)

      integer :: i, l0, l1, l2, l3, l4, tile_count, u_size2
      real*8 :: error, f0, f1, u_old

      u_size2 = u_size*u_size

      error = eps + 1.0d0
      tile_count = (u_size - 2)/tile_size
      f0 = 1.0d0 - omega
      f1 = omega / 6.0d0

      iter = 0

      do while (error > eps .and. iter < max_iter)

         i = u_size*(u_size + 1) + 2
         error = 0.0d0

         do l4 = 3, u_size
            do l3 = 1, tile_count
               do l2 = 1, tile_count
                  do l1 = 1, tile_size
                     do l0 = 1, tile_size

                        u_old = u(i)
                        u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                           + u(i+u_size2) + u(i-u_size2))

                        error = error + abs(u(i) - u_old)
                        i = i + 1

                     end do
                     i = i + u_size - tile_size
                  end do
                  i = i - tile_size*(u_size - 1)
               end do
               i = i + u_size*(tile_size - 1) + 2
            end do
            i = i + u_size*2
         end do
         iter = iter + 1
      end do

   end subroutine tiling_3d_sor

   !SOR с тестовым субтайлингом.
   subroutine subtiling_3d_sor_test(u_1d, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout), target :: u_1d(u_size*u_size*u_size)
      real*8, pointer :: u(:,:,:)
      integer :: l0, l1, l2, l3, l4, l5, l6
      real*8 :: f0, f1, u_old, error

      u(1:u_size, 1:u_size, 1:u_size) => u_1d
      
      f0 = 1.0d0 - omega
      f1 = omega / 6.0d0

      do iter = 1, max_iter

         error = 0.0d0

         do l6 = 2, u_size, tile_size
            do l5 = 2, u_size, tile_size
               do l4 = 2, u_size, tile_size
                  do l3 = 0, tile_size - 1
                     do l2 = max(l6 - l3,2), min(l6 + tile_size - l3 - 1,u_size-1)
                        do l1 = max(l5 - l3,2), min(l5 + tile_size - l3 - 1,u_size-1)
                           do l0 = max(l4 - l3,2), min(l4 + tile_size - l3 - 1,u_size-1)

                              u_old = u(l0,l1,l2)
                              u(l0,l1,l2) = f0*u_old + f1*(u(l0-1,l1,l2) + u(l0+1,l1,l2) + &
                                 u(l0,l1-1,l2) + u(l0,l1+1,l2) + u(l0,l1,l2-1) + u(l0,l1,l2+1))
                              if (l3.eq.tile_size-1) then
                                 error = error + abs(u(l0,l1,l2) - u_old)
                              end if

                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

         if (error < eps) then
            return
         end if

      end do
   end subroutine subtiling_3d_sor_test

   !SOR с субтайлингом.
   subroutine subtiling_3d_sor(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8, intent(inout) :: u(u_size*u_size*u_size)

      integer :: tile_count, i, l0, l1, l2, l3, l4, l5, u_size2
      real*8 :: f0, f1, u_old, error

      u_size2 = u_size*u_size

      tile_count = (u_size - 2)/tile_size

      f0 = 1.0d0 - omega
      f1 = omega / 6.0d0

      do iter = 1, max_iter

         error = 0.0d0
         i = u_size*(u_size + 1) + 2

         do l5 = 3, u_size
            do l2 = 0, subtile_level
               do l1 = 1, tile_size - l2
                  do l0 = 1, tile_size - l2

                     u_old = u(i)

                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                        + u(i+u_size2) + u(i+u_size2))

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
                        u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                           + u(i+u_size2) + u(i+u_size2))

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
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                        + u(i+u_size2) + u(i+u_size2))

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
                        u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                           + u(i+u_size2) + u(i+u_size2))

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
                           u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                              + u(i+u_size2) + u(i+u_size2))

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
                        u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                           + u(i+u_size2) + u(i+u_size2))

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
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                        + u(i+u_size2) + u(i+u_size2))

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
                        u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                           + u(i+u_size2) + u(i+u_size2))

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
                     u(i) = f0*u_old + f1*(u(i-1) + u(i+1) + u(i-u_size) + u(i+u_size) &
                        + u(i+u_size2) + u(i+u_size2))

                     if (l2.eq.subtile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(l1+1) + l0 + 2
         end do

         if (error < eps) then
            return
         end if

      end do
   end subroutine subtiling_3d_sor
end module subgrid_3d_solver