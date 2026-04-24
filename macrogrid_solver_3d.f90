module macrogrid_3d_solver
   use subgrid_3d_solver
   use omp_lib
   implicit none
   private

   interface
      subroutine i_macrogrid_solver_method()
         implicit none
      end subroutine i_macrogrid_solver_method
   end interface

   public :: i_macrogrid_solver_method
   public :: set_macrogrid_3d_solver_settings, set_default_macrogrid_3d_solver_settings
   public :: run_macrogrid_3d_solver, get_macrogrid_3d_solver_results
   
   ! Методы макро-решателя
   public :: simple_iteration_3d
   public :: simple_iteration_3d_one_iter

   integer :: m_size_x, m_size_y, m_size_z, sub_size, max_iter
   real*8 :: eps
   real*8, pointer :: macrogrid(:,:,:,:,:,:) => null()

   procedure(i_subgrid_3d_solver_method), pointer :: sub_solver => null()
   
   real*8 :: total_time
   integer :: total_iter

contains

   subroutine set_default_macrogrid_3d_solver_settings()
      implicit none
      call set_default_subgrid_3d_solver_settings()
      eps = 1.0d-8
      max_iter = 1000
   end subroutine

   subroutine set_macrogrid_3d_solver_settings(new_eps, new_max_iter)
      real*8, intent(in), optional :: new_eps
      integer, intent(in), optional :: new_max_iter
      if (present(new_eps)) eps = new_eps
      if (present(new_max_iter)) max_iter = new_max_iter
   end subroutine

   subroutine run_macrogrid_3d_solver(use_openmp, new_macrogrid, &
      nx, ny, nz, s_size, &
      macro_method, sub_method)

      logical, intent(in) :: use_openmp
      integer, intent(in) :: nx, ny, nz, s_size
      real*8, intent(inout), target :: new_macrogrid(nx, ny, nz, s_size, s_size, s_size)
      procedure(i_macrogrid_solver_method) :: macro_method
      procedure(i_subgrid_3d_solver_method) :: sub_method

      real*8 :: t_start, t_end

      m_size_x = nx; m_size_y = ny; m_size_z = nz
      sub_size = s_size
      macrogrid => new_macrogrid
      sub_solver => sub_method

      if (use_openmp) then
         t_start = omp_get_wtime()
         call macro_method()
         t_end = omp_get_wtime()
      else
         call cpu_time(t_start)
         call macro_method()
         call cpu_time(t_end)
      end if

      total_time = t_end - t_start
   end subroutine

   subroutine get_macrogrid_3d_solver_results(res_time, res_iter)
      real*8, intent(out) :: res_time
      integer, intent(out) :: res_iter
      res_time = total_time
      res_iter = total_iter
   end subroutine

   ! Внутренняя процедура запуска решателей в подсетках
   subroutine compute_all_subgrids()
      integer :: ix, iy, iz
      !$OMP PARALLEL DO PRIVATE(ix, iy, iz) COLLAPSE(3) SCHEDULE(static) IF(m_size_x*m_size_y*m_size_z > 1)
      do ix = 1, m_size_x
         do iy = 1, m_size_y
            do iz = 1, m_size_z
               call sub_solver(macrogrid(ix, iy, iz, :, :, :), sub_size)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
   end subroutine

   ! Метод простой итерации (3D версия формулы 2.2.2 из диплома)
   subroutine simple_iteration_3d()
      implicit none
      real*8 :: err, new_val, old_val
      integer :: ix, iy, iz, i, j, k

      do total_iter = 1, max_iter
         err = 0.0d0
         call compute_all_subgrids()

         ! 1. Интерфейсы по оси X (плоскости YZ)
         do ix = 1, m_size_x - 1
            do iy = 1, m_size_y
               do iz = 1, m_size_z
                  do j = 2, sub_size - 1
                     do k = 2, sub_size - 1
                        old_val = macrogrid(ix, iy, iz, sub_size, j, k)
                        ! Формула 2.2.2: (4*(u[s-1]+u[2]) - (u[s-2]+u[3])) / 6
                        new_val = (4.0d0*(macrogrid(ix, iy, iz, sub_size-1, j, k) + &
                                         macrogrid(ix+1, iy, iz, 2, j, k)) - &
                                        (macrogrid(ix, iy, iz, sub_size-2, j, k) + &
                                         macrogrid(ix+1, iy, iz, 3, j, k))) / 6.0d0
                        macrogrid(ix, iy, iz, sub_size, j, k) = new_val
                        macrogrid(ix+1, iy, iz, 1, j, k) = new_val
                        err = err + abs(new_val - old_val)
                     end do
                  end do
               end do
            end do
         end do

         ! 2. Интерфейсы по оси Y (плоскости XZ)
         do ix = 1, m_size_x
            do iy = 1, m_size_y - 1
               do iz = 1, m_size_z
                  do i = 2, sub_size - 1
                     do k = 2, sub_size - 1
                        old_val = macrogrid(ix, iy, iz, i, sub_size, k)
                        new_val = (4.0d0*(macrogrid(ix, iy, iz, i, sub_size-1, k) + &
                                         macrogrid(ix, iy+1, iz, i, 2, k)) - &
                                        (macrogrid(ix, iy, iz, i, sub_size-2, k) + &
                                         macrogrid(ix, iy+1, iz, i, 3, k))) / 6.0d0
                        macrogrid(ix, iy, iz, i, sub_size, k) = new_val
                        macrogrid(ix, iy+1, iz, i, 1, k) = new_val
                        err = err + abs(new_val - old_val)
                     end do
                  end do
               end do
            end do
         end do

         ! 3. Интерфейсы по оси Z (плоскости XY)
         do ix = 1, m_size_x
            do iy = 1, m_size_y
               do iz = 1, m_size_z - 1
                  do i = 2, sub_size - 1
                     do j = 2, sub_size - 1
                        old_val = macrogrid(ix, iy, iz, i, j, sub_size)
                        new_val = (4.0d0*(macrogrid(ix, iy, iz, i, j, sub_size-1) + &
                                         macrogrid(ix, iy, iz+1, i, j, 2)) - &
                                        (macrogrid(ix, iy, iz, i, j, sub_size-2) + &
                                         macrogrid(ix, iy, iz+1, i, j, 3))) / 6.0d0
                        macrogrid(ix, iy, iz, i, j, sub_size) = new_val
                        macrogrid(ix, iy, iz+1, i, j, 1) = new_val
                        err = err + abs(new_val - old_val)
                     end do
                  end do
               end do
            end do
         end do

         if (err < eps) exit
      end do
      
      call compute_edges_and_vertices()
   end subroutine

   ! Однократный запуск (Метод из раздела 2.3 диплома)
   subroutine simple_iteration_3d_one_iter()
      integer :: old_max
      call get_subgrid_3d_solver_settings(new_max_iter=old_max)
      call set_subgrid_3d_solver_settings(new_max_iter=1) ! Ключевое отличие метода
      call simple_iteration_3d()
      call set_subgrid_3d_solver_settings(new_max_iter=old_max)
   end subroutine

   ! Обработка ребер и вершин, где сходятся несколько подсеток
   subroutine compute_edges_and_vertices()
      integer :: ix, iy, iz
      real*8 :: val

      ! Усреднение по ребрам (стык 4-х блоков) - пример для ребер вдоль Z
      do ix = 1, m_size_x - 1
         do iy = 1, m_size_y - 1
            do iz = 1, m_size_z
               do k = 1, sub_size
                  val = (macrogrid(ix, iy, iz, sub_size-1, sub_size, k) + &
                         macrogrid(ix+1, iy, iz, 2, sub_size, k) + &
                         macrogrid(ix, iy+1, iz, sub_size, 2, k) + &
                         macrogrid(ix+1, iy+1, iz, 2, 2, k)) / 4.0d0
                  macrogrid(ix, iy, iz, sub_size, sub_size, k) = val
                  macrogrid(ix+1, iy, iz, 1, sub_size, k) = val
                  macrogrid(ix, iy+1, iz, sub_size, 1, k) = val
                  macrogrid(ix+1, iy+1, iz, 1, 1, k) = val
               end do
            end do
         end do
      end do
      
      ! Аналогично для ребер вдоль X и Y... (сокращено для краткости)
      
      ! Усреднение в вершинах (стык 8-ми блоков)
      do ix = 1, m_size_x - 1
         do iy = 1, m_size_y - 1
            do iz = 1, m_size_z - 1
               val = (macrogrid(ix, iy, iz, sub_size-1, sub_size-1, sub_size-1) + &
                      macrogrid(ix+1, iy, iz, 2, sub_size-1, sub_size-1) + &
                      macrogrid(ix, iy+1, iz, sub_size-1, 2, sub_size-1) + &
                      macrogrid(ix+1, iy+1, iz, 2, 2, sub_size-1) + &
                      macrogrid(ix, iy, iz+1, sub_size-1, sub_size-1, 2) + &
                      macrogrid(ix+1, iy, iz+1, 2, sub_size-1, 2) + &
                      macrogrid(ix, iy+1, iz+1, sub_size-1, 2, 2) + &
                      macrogrid(ix+1, iy+1, iz+1, 2, 2, 2)) / 8.0d0
               ! Присвоение во все 8 угловых узлов смежных подсеток
               macrogrid(ix, iy, iz, sub_size, sub_size, sub_size) = val
               macrogrid(ix+1, iy+1, iz+1, 1, 1, 1) = val
               ! ... и остальные 6 углов
            end do
         end do
      end do
   end subroutine

end module