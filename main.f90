program main
   use macrogrid_boundary_initializer
   use subgrid_boundary_initializer
   use subgrid_3d_boundary_initializer
   use macrogrid_solver
   use subgrid_solver
   use subgrid_3d_solver

   implicit none
   ! io: номер файлового дескриптора для записи результатов экспериментов.
   integer, parameter :: io = 10

   ! subgrid_sizes: набор тестовых размеров двумерных и трехмерных подсеток.
   integer, dimension(8), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]

   ! tile_sizes: набор размеров плиток для tiled- и subtiling-вариантов.
   integer, dimension(4), parameter :: tile_sizes = &
      [1, 2, 4, 8]

   ! subgrid_omegas: подобранные значения omega для двумерного оригинального SOR.
   real*8, dimension(8), parameter :: subgrid_omegas = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]
   ! subgrid_const_omegas: подобранные omega для двумерной задачи с постоянной границей.
   real*8, dimension(8), parameter :: subgrid_const_omegas = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]
   ! subgrid_3d_const_omegas: подобранные omega для трехмерной задачи с постоянной границей.
   real*8, dimension(8), parameter :: subgrid_3d_const_omegas = &
      [1.49519d0, 1.69289d0, 1.82788d0, 1.90853d0, 1.95281d0, 1.97609d0, 1.5d0, 1.5d0]

   ! subgrid_repeats_count: число повторов двумерных запусков для усреднения времени.
   integer, dimension(8), parameter :: subgrid_repeats_count = &
      [10000, 10000, 1000, 1000, 100, 100, 10, 10]
   ! subgrid_3d_repeats_count: число повторов трехмерных запусков для усреднения времени.
   integer, dimension(8), parameter :: subgrid_3d_repeats_count = &
      [10000, 10000, 1000, 100, 10, 1, 1, 1]

   ! macrogrid_sizes: набор размеров макросетки по каждой координате.
   integer, dimension(8), parameter :: macrogrid_sizes = &
      [1, 2, 3, 4, 5, 6, 7, 8]
   ! macrogrid_omegas_1, macrogrid_omegas_2, macrogrid_omegas_4: заранее подобранные omega
   ! для разных размеров макросетки и размеров подсетки.
   real*8, dimension(8), parameter :: macrogrid_omegas_1 = &
      [1.9d0, 1.99d0, 1.999d0, 1.999d0, 1.999d0, 1.999d0, 1.999d0, 1.999d0]
   real*8, dimension(8), parameter :: macrogrid_omegas_2 = &
      [1.97137756696708d0, 1.98302605323770d0, 1.99130449164847d0, &
      1.99618167072694d0, 1.99801933675286d0, 1.99950121735996d0, &
      1.99926204012712d0, 1.999d0]
   real*8, dimension(8), parameter :: macrogrid_omegas_4 = &
      [1.99283687348090d0, 1.99684870120073d0, 1.99836547729607d0, &
      1.99866111661454d0, 1.99955171858458d0, 1.999d0, &
      1.999d0, 1.999d0]
   ! macrogrid_omegas: сводная таблица omega для макросеточных экспериментов.
   real*8, dimension(8,8) :: macrogrid_omegas

   macrogrid_omegas(1,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(2,:) = macrogrid_omegas_2(:)
   macrogrid_omegas(3,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(4,:) = macrogrid_omegas_4(:)
   macrogrid_omegas(5,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(6,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(7,:) = macrogrid_omegas_1(:)
   macrogrid_omegas(8,:) = macrogrid_omegas_1(:)

   call subgrid_3d_test()

contains
   ! Запускает серию тестов для трехмерных решателей подсетки.
   subroutine subgrid_3d_test()
      implicit none

      integer, dimension(6), parameter :: subgrid_configs = &
         [1, 2, 3, 4, 5, 6]

      integer, dimension(3), parameter :: tile_configs = &
         [2, 3, 4]

      integer :: sub_cfg, tile_cfg

      write(*, *) "Running subgrid 3D tests"

      open(io, file='results.txt', status='unknown', action='write', recl=10000)

      write(io, *) "original_3d_sor"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_3d_test(subgrid_solver_method = original_3d_sor, &
            subgrid_size_idx = subgrid_configs(sub_cfg), tile_size_idx = 1)
      end do
      write(io, *) " "

      write(io, *) "tiling_3d_sor"
      do tile_cfg = 1, size(tile_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_subgrid_3d_test(subgrid_solver_method = tiling_3d_sor, &
               subgrid_size_idx = subgrid_configs(sub_cfg), tile_size_idx = tile_configs(tile_cfg))
         end do
      end do
      write(io, *) " "

      ! write(io, *) "subtiling_3d_sor"
      ! do tile_cfg = 1, size(tile_configs)
      !    do sub_cfg = 1, size(subgrid_configs)
      !       call run_subgrid_3d_test(subgrid_solver_method = subtiling_3d_sor_test, &
      !          subgrid_size_idx = subgrid_configs(sub_cfg), tile_size_idx = tile_configs(tile_cfg))
      !    end do
      ! end do
      ! write(io, *) " "

      close(io)

      write(*, *) "End of subgrid tests"
   end subroutine subgrid_3d_test
   ! Запускает серию тестов для двумерных решателей подсетки.
   subroutine subgrid_test()
      implicit none

      integer, dimension(8), parameter :: subgrid_configs = &
         [1, 2, 3, 4, 5, 6, 7, 8]

      integer :: sub_cfg

      write(*, *) "Running subgrid tests"

      open(io, file='results.txt', status='unknown', action='write', recl=10000)

      write(io, *) "original_sor"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = original_sor, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "tiling_sor_1"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = tiling_sor_1, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "tiling_sor_2"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = tiling_sor_2, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "tiling_sor_4"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = tiling_sor_4, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "tiling_sor_8"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = tiling_sor_8, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "tiling_sor_16"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = tiling_sor_16, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_1"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = subtiling_sor_1, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_2"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = subtiling_sor_2, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_4"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = subtiling_sor_4, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_8"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = subtiling_sor_8, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      write(io, *) "subtiling_sor_16"
      do sub_cfg = 1, size(subgrid_configs)
         call run_subgrid_test(subgrid_solver_method = subtiling_sor_16, &
            subgrid_size_idx = subgrid_configs(sub_cfg))
      end do
      write(io, *) " "

      close(io)

      write(*, *) "End of subgrid tests"
   end subroutine subgrid_test
   ! Запускает серию тестов для макросеточного решателя.
   subroutine macrogrid_test()
      implicit none

      integer, dimension(4), parameter :: subgrid_configs = &
         [1, 2, 3, 4]

      integer, dimension(1), parameter :: macrogrid_configs = &
         [2]

      integer :: mac_cfg, sub_cfg

      write(*, *) "Running macrogrid tests"

      open(io, file='results.txt', status='unknown', action='write', recl=10000)

      write(io, *) "sor_fixed_omega"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .false., &
               macrogrid_solver_method = sor_fixed_omega, subgrid_solver_method = original_sor, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "sor_fixed_omega & one_iter"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .false., &
               macrogrid_solver_method = sor_fixed_omega_one_iter, subgrid_solver_method = original_sor, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "sor_fixed_omega & one_iter & tiling_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .false., &
               macrogrid_solver_method = sor_fixed_omega_one_iter, subgrid_solver_method = tiling_sor_8, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "sor_fixed_omega & one_iter & tiling_sor & openmp"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .true., &
               macrogrid_solver_method = sor_fixed_omega_one_iter, subgrid_solver_method = tiling_sor_8, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .false., &
               macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = original_sor, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals & tiling_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .false., &
               macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = tiling_sor_8, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals & subtiling_sor"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .false., &
               macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = subtiling_sor_4, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      write(io, *) "conjugate_residuals & subtiling_sor & openmp"
      do mac_cfg = 1, size(macrogrid_configs)
         do sub_cfg = 1, size(subgrid_configs)
            call run_macrogrid_test(use_openmp = .true., &
               macrogrid_solver_method = conjugate_residuals, subgrid_solver_method = subtiling_sor_4, &
               macrogrid_size_idx = macrogrid_configs(mac_cfg), subgrid_size_idx = subgrid_configs(sub_cfg))
         end do
      end do
      write(io, *) " "

      close(io)

      write(*, *) "End of macrogrid tests"
   end subroutine macrogrid_test
   ! Выполняет один набор измерений для выбранного трехмерного решателя подсетки.
   ! Аргументы:
   ! subgrid_solver_method: тестируемый трехмерный алгоритм решения подсетки.
   ! subgrid_size_idx: индекс размера подсетки в таблице subgrid_sizes.
   ! tile_size_idx: индекс размера плитки в таблице tile_sizes.
   subroutine run_subgrid_3d_test(subgrid_solver_method, subgrid_size_idx, tile_size_idx)
      implicit none
      procedure(i_subgrid_3d_solver_method) :: subgrid_solver_method
      integer, intent(in) :: subgrid_size_idx, tile_size_idx

      real*8, allocatable :: subgrid(:)

      real*8 :: error, time, sum_time
      integer :: iter, i

      call set_default_subgrid_3d_solver_settings()
      call set_subgrid_3d_solver_settings(new_omega = subgrid_3d_const_omegas(subgrid_size_idx), &
         new_eps = 1.0d-5, new_tile_size = tile_sizes(tile_size_idx), &
         new_subtile_level = tile_sizes(tile_size_idx))

      allocate(subgrid(subgrid_sizes(subgrid_size_idx)&
         *subgrid_sizes(subgrid_size_idx)*subgrid_sizes(subgrid_size_idx)))

      sum_time = 0.0d0

      do i = 1, subgrid_3d_repeats_count(subgrid_size_idx)

         call initialize_subgrid_3d_boundary(subgrid, subgrid_sizes(subgrid_size_idx))

         call run_subgrid_3d_solver(subgrid, subgrid_sizes(subgrid_size_idx), subgrid_solver_method)

         call get_subgrid_3d_solver_results(time, iter)

         sum_time = sum_time + time

      end do

      call compute_subgrid_3d_boundary_error(subgrid, subgrid_sizes(subgrid_size_idx), error)

      time = sum_time / real(subgrid_3d_repeats_count(subgrid_size_idx))

      write(io,*) "Tile_size&Subgrid_size&Error&Run_time&Iterations#", &
         tile_sizes(tile_size_idx), subgrid_sizes(subgrid_size_idx), "&", error, "&", time, "&", iter

      deallocate(subgrid)

   end subroutine run_subgrid_3d_test
   ! Выполняет один набор измерений для выбранного двумерного решателя подсетки.
   ! Аргументы:
   ! subgrid_solver_method: тестируемый двумерный алгоритм решения подсетки.
   ! subgrid_size_idx: индекс размера подсетки в таблице subgrid_sizes.
   subroutine run_subgrid_test(subgrid_solver_method, subgrid_size_idx)
      implicit none
      procedure(i_subgrid_solver_method) :: subgrid_solver_method
      integer, intent(in) :: subgrid_size_idx

      real*8, allocatable :: subgrid(:)

      real*8 :: error, time, sum_time
      integer :: iter, i

      call set_default_subgrid_solver_settings()
      call set_subgrid_solver_settings(new_omega = subgrid_omegas(subgrid_size_idx), new_eps = 1.0d-5)

      allocate(subgrid(subgrid_sizes(subgrid_size_idx)*subgrid_sizes(subgrid_size_idx)))

      sum_time = 0.0d0

      do i = 1, subgrid_repeats_count(subgrid_size_idx)

         call initialize_subgrid_boundary(subgrid, subgrid_sizes(subgrid_size_idx))

         call run_subgrid_solver(subgrid, subgrid_sizes(subgrid_size_idx), subgrid_solver_method)

         call get_subgrid_solver_results(time, iter)

         sum_time = sum_time + time

      end do

      call compute_subgrid_boundary_error(subgrid, subgrid_sizes(subgrid_size_idx), error)

      time = sum_time / real(subgrid_repeats_count(subgrid_size_idx))

      write(io,*) "Subgrid_size&Error&Run_time&Iterations#", &
         subgrid_sizes(subgrid_size_idx), "&", error, "&", time, "&", iter

      deallocate(subgrid)

   end subroutine run_subgrid_test
   ! Выполняет один набор измерений для выбранной комбинации макро- и подрешателей.
   ! Аргументы:
   ! use_openmp: признак использования OpenMP для расчета подсеток.
   ! macrogrid_solver_method: тестируемый алгоритм согласования интерфейсов.
   ! subgrid_solver_method: тестируемый решатель отдельной подсетки.
   ! macrogrid_size_idx: индекс размера макросетки в таблице macrogrid_sizes.
   ! subgrid_size_idx: индекс размера подсетки в таблице subgrid_sizes.
   subroutine run_macrogrid_test(use_openmp, macrogrid_solver_method, subgrid_solver_method, &
      macrogrid_size_idx, subgrid_size_idx)
      implicit none
      procedure(i_macrogrid_solver_method) :: macrogrid_solver_method
      procedure(i_subgrid_solver_method) :: subgrid_solver_method
      logical, intent(in) :: use_openmp
      integer, intent(in) :: macrogrid_size_idx, subgrid_size_idx

      real*8, allocatable :: macrogrid(:, :, :, :)

      real*8 :: error, time
      integer :: iter

      call set_default_macrogrid_solver_settings()
      call set_subgrid_solver_settings(new_omega = subgrid_omegas(subgrid_size_idx))
      call set_macrogrid_solver_settings(new_omega = macrogrid_omegas(macrogrid_size_idx, subgrid_size_idx))

      allocate(macrogrid(macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), subgrid_sizes(subgrid_size_idx)))

      call initialize_macrogrid_boundary(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))

      call run_macrogrid_solver(use_openmp, macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), macrogrid_solver_method, subgrid_solver_method)

      call get_macrogrid_solver_results(time, iter)

      call compute_macrogrid_boundary_error(macrogrid, &
         macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), error)

      write(io,*) "Macrogrid_size&Subgrid_size&Error&Run_time&Iterations#", &
         macrogrid_sizes(macrogrid_size_idx), "&", subgrid_sizes(subgrid_size_idx), &
         "&", error, "&", time, "&", iter

      deallocate(macrogrid)

   end subroutine run_macrogrid_test
   ! Подбирает omega для трехмерного решателя методом золотого сечения.
   ! Аргументы:
   ! subgrid_size_idx: индекс размера подсетки в таблице subgrid_sizes.
   ! subgrid_omega: найденное оптимальное значение omega.
   subroutine subgrid_3d_golden_section(subgrid_size_idx, subgrid_omega)
      implicit none
      integer, intent(in) :: subgrid_size_idx
      real*8, intent(out) :: subgrid_omega

      real*8 :: eps_subgrid = 1.0d-5, eps_golden = 1.0d-5
      integer :: iter_interface = 10000

      real*8, allocatable :: subgrid(:)

      real*8 :: error, time, time1, time2
      integer :: iter,iter1,iter2

      real*8 :: k, a, b, wa, wb

      call set_default_subgrid_3d_solver_settings()
      call set_subgrid_3d_solver_settings(new_omega = subgrid_3d_const_omegas(subgrid_size_idx), &
         new_eps = 1.0d-5, new_tile_size = 4, new_subtile_level = 4, new_max_iter = iter_interface)

      allocate(subgrid(subgrid_sizes(subgrid_size_idx)&
         *subgrid_sizes(subgrid_size_idx)*subgrid_sizes(subgrid_size_idx)))

      k = 0.6180339887d0
      a = 1.0d0
      b = 2.0d0
      wa = 1.0d0
      wb = 2.0d0

      subgrid_omega = 1.5d0

      do
         if (abs(wa-wb) < eps_golden) then
            subgrid_omega = (wb+wa)/2
            exit
         end if

         wa = a + (1-k)*(b-a)
         wb = a + k*(b-a)

         call initialize_subgrid_3d_boundary(subgrid, subgrid_sizes(subgrid_size_idx))
         call set_subgrid_3d_solver_settings(new_omega = wa)
         call run_subgrid_3d_solver(subgrid, subgrid_sizes(subgrid_size_idx), tiling_3d_sor)
         call get_subgrid_3d_solver_results(time, iter)

         iter1 = iter
         time1 = time
         time = 0.0d0
         iter = 0

         call initialize_subgrid_3d_boundary(subgrid, subgrid_sizes(subgrid_size_idx))
         call set_subgrid_3d_solver_settings(new_omega = wb)
         call run_subgrid_3d_solver(subgrid, subgrid_sizes(subgrid_size_idx), tiling_3d_sor)
         call get_subgrid_3d_solver_results(time, iter)

         iter2 = iter
         time2 = time
         time = 0.0d0
         iter = 0

         if (iter1 > iter2) then
            a = wa
         else
            b = wb
         end if

         write(io,*) "a = ", a, "b = ", b

      end do

      deallocate(subgrid)

   end subroutine subgrid_3d_golden_section
   ! Подбирает omega для макросеточного решателя методом золотого сечения.
   ! Аргументы:
   ! macrogrid_size_idx: индекс размера макросетки в таблице macrogrid_sizes.
   ! subgrid_size_idx: индекс размера подсетки в таблице subgrid_sizes.
   ! macrogrid_omega: найденное оптимальное значение omega.
   subroutine golden_section(macrogrid_size_idx, subgrid_size_idx, macrogrid_omega)
      implicit none
      integer, intent(in) :: macrogrid_size_idx, subgrid_size_idx
      real*8, intent(out) :: macrogrid_omega

      real*8 :: eps_subgrid = 1.0d-8, eps_interface = 1.0d-8, eps_golden = 1.0d-5
      integer :: max_iter_interface = 100000, max_iter_subgrid = 1

      real*8, allocatable :: macrogrid(:, :, :, :)

      real*8 :: error, time, time1, time2
      integer :: iter,iter1,iter2

      real*8 :: k, a, b, wa, wb

      call set_default_macrogrid_solver_settings()
      call set_subgrid_solver_settings(new_omega = subgrid_omegas(subgrid_size_idx))

      k = 0.6180339887d0
      a = 1.0d0
      b = 2.0d0
      wa = 1.0d0
      wb = 2.0d0

      allocate(macrogrid(macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
         subgrid_sizes(subgrid_size_idx), subgrid_sizes(subgrid_size_idx)))

      macrogrid_omega = 1.5d0

      do
         if (abs(wa-wb) < eps_golden) then
            macrogrid_omega = (wb+wa)/2
            exit
         end if

         wa = a + (1-k)*(b-a)
         wb = a + k*(b-a)

         call initialize_macrogrid_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
         call set_macrogrid_solver_settings(new_omega = wa)
         call run_macrogrid_solver(.false., macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
         call get_macrogrid_solver_results(time, iter)
         iter1 = iter
         time1 = time
         time = 0.0d0
         iter = 0

         call initialize_macrogrid_boundary(macrogrid, &
            macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), subgrid_sizes(subgrid_size_idx))
         call set_macrogrid_solver_settings(new_omega = wb)
         call run_macrogrid_solver(.false., macrogrid, macrogrid_sizes(macrogrid_size_idx), macrogrid_sizes(macrogrid_size_idx), &
            subgrid_sizes(subgrid_size_idx), sor_fixed_omega, tiling_sor)
         call get_macrogrid_solver_results(time, iter)
         iter2 = iter
         time2 = time
         time = 0.0d0
         iter = 0

         if (iter1 > iter2) then
            a = wa
         else
            b = wb
         end if

      end do

      deallocate(macrogrid)

   end subroutine golden_section
   ! Печатает значения макросетки в удобном для просмотра виде.
   ! Аргументы:
   ! macrogrid: четырехмерный массив значений макросетки.
   ! macrogrid_size: число подсеток по каждой координате.
   ! subgrid_size: размер одной подсетки.
   subroutine print_macrogrid(macrogrid, macrogrid_size, subgrid_size)
      implicit none
      integer, intent(in) :: macrogrid_size, subgrid_size
      real*8, intent(in) :: macrogrid(macrogrid_size, macrogrid_size, subgrid_size, subgrid_size)
      integer :: i, j, x, y

      do i = 1, macrogrid_size
         do y = 1, subgrid_size
            do j = 1, macrogrid_size
               do x = 1, subgrid_size
                  write(*, '(10F8.3)', advance='no') macrogrid(j, i, x, y)
               end do
            end do
            print *, ' '
         end do
      end do
      print *, ' '

   end subroutine print_macrogrid
   ! Запускает подбор omega для набора трехмерных размеров подсетки.
   subroutine subgrid_3d_golden_section_test()
      implicit none

      integer, dimension(8), parameter :: subgrid_configs = &
         [1, 2, 3, 4, 5, 6, 7, 8]

      integer :: sub_cfg
      real*8 :: res

      write(*, *) "Running subgrid 3D golden section tests"

      open(io, file='results.txt', status='unknown', action='write', recl=10000)

      do sub_cfg = 1, size(subgrid_configs)
         write(io, *) "Subgrid 3D golden section results:"
         write(io, *) "subgrid_size_idx = ", subgrid_configs(sub_cfg)
         call subgrid_3d_golden_section(subgrid_size_idx = subgrid_configs(sub_cfg), subgrid_omega = res)
         write(io, *) "Omega = ", res
         write(io, *) " "
      end do

      close(io)

      write(*, *) "End of subgrid tests"
   end subroutine subgrid_3d_golden_section_test
end program main
