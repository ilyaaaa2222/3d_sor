module macrogrid_solver
   use subgrid_solver
   use omp_lib
   implicit none
   private

   interface
      ! Интерфейс макроитерационного метода, который работает с уже настроенным модулем.
      subroutine i_macrogrid_solver_method()
         implicit none
      end subroutine i_macrogrid_solver_method

      ! Интерфейс внутренней процедуры вычисления всех подсеток макросетки.
      subroutine i_subgrids_computer()
         implicit none
      end subroutine i_subgrids_computer
   end interface

   public :: i_macrogrid_solver_method, i_subgrid_solver_method
   public :: set_macrogrid_solver_settings, set_default_macrogrid_solver_settings
   public :: get_macrogrid_solver_settings, get_macrogrid_solver_results
   public :: run_macrogrid_solver

   public :: simple_iteration
   public :: simple_iteration_one_iter
   public :: sor
   public :: sor_fixed_omega
   public :: sor_fixed_omega_one_iter
   public :: conjugate_residuals

   ! Параметры и состояние макрорешателя.
   ! macrogrid_size_x, macrogrid_size_y: размеры разбиения области на подсетки по осям X и Y.
   ! subgrid_size: размер одной квадратной подсетки по одной координате.
   ! max_iter: верхний предел числа макроитераций.
   integer :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter
   ! omega: параметр релаксации для SOR-методов на интерфейсах.
   ! eps: критерий остановки по суммарной невязке интерфейсов.
   real*8 :: omega, eps
   ! macrogrid: указатель на рабочий массив значений всей макросетки.
   real*8, pointer :: macrogrid(:,:,:,:)

   ! dx, dy: глобальные шаги сетки по осям X и Y.
   real*8 :: dx, dy
   ! interface_size: число внутренних интерфейсных узлов между подсетками.
   integer :: interface_size
   ! subgrid_solver_method: выбранный решатель отдельной подсетки.
   procedure(i_subgrid_solver_method), pointer :: subgrid_solver_method => null()
   ! subgrids_computer: выбранный способ обхода подсеток, последовательный или OpenMP.
   procedure(i_subgrids_computer), pointer :: subgrids_computer => null()

   ! time: время последнего запуска макрорешателя.
   real*8 :: time
   ! iter: число итераций, выполненных последним запуском.
   integer :: iter

contains
   ! Устанавливает значения по умолчанию для макрорешателя и вложенного решателя подсеток.
   subroutine set_default_macrogrid_solver_settings()
      implicit none

      call set_default_subgrid_solver_settings()

      eps = 1.0d-8
      max_iter = 100000
      omega = -1

   end subroutine set_default_macrogrid_solver_settings
   ! Обновляет параметры макрорешателя.
   ! Аргументы:
   ! new_eps: новое значение критерия остановки.
   ! new_max_iter: новое ограничение на число макроитераций.
   ! new_omega: новое значение параметра релаксации.
   subroutine set_macrogrid_solver_settings(new_eps, new_max_iter, new_omega)
      implicit none
      real*8, intent(in), optional :: new_eps, new_omega
      integer, intent(in), optional :: new_max_iter

      if (present(new_eps)) eps = new_eps
      if (present(new_max_iter)) max_iter = new_max_iter
      if (present(new_omega)) omega = new_omega

   end subroutine set_macrogrid_solver_settings
   ! Возвращает текущие настройки макрорешателя.
   ! Аргументы:
   ! new_eps: сюда записывается текущее значение критерия остановки.
   ! new_max_iter: сюда записывается текущее ограничение на число итераций.
   ! new_omega: сюда записывается текущее значение параметра релаксации.
   subroutine get_macrogrid_solver_settings(new_eps, new_max_iter, new_omega)
      implicit none
      real*8, intent(out), optional :: new_eps, new_omega
      integer, intent(out), optional :: new_max_iter

      if (present(new_eps)) new_eps = eps
      if (present(new_max_iter)) new_max_iter = max_iter
      if (present(new_omega)) new_omega = omega

   end subroutine get_macrogrid_solver_settings
   ! Настраивает модуль, запускает выбранный макрометод и измеряет время работы.
   ! Аргументы:
   ! use_openmp: признак параллельного вычисления подсеток через OpenMP.
   ! new_macrogrid: рабочий массив макросетки, который будет обновляться in-place.
   ! new_macrogrid_size_x: число подсеток по оси X.
   ! new_macrogrid_size_y: число подсеток по оси Y.
   ! new_subgrid_size: размер одной квадратной подсетки.
   ! new_macrogrid_solver_method: выбранный метод согласования интерфейсов между подсетками.
   ! new_subgrid_solver_method: выбранный метод решения отдельной подсетки.
   subroutine run_macrogrid_solver(use_openmp, new_macrogrid, &
      new_macrogrid_size_x, new_macrogrid_size_y, new_subgrid_size, &
      new_macrogrid_solver_method, new_subgrid_solver_method)

      implicit none
      logical, intent(in) :: use_openmp
      integer, intent(in) :: new_macrogrid_size_x, new_macrogrid_size_y, new_subgrid_size
      real*8, intent(inout), target :: new_macrogrid(:,:,:,:)
      procedure(i_macrogrid_solver_method) :: new_macrogrid_solver_method
      procedure(i_subgrid_solver_method) :: new_subgrid_solver_method

      real*8 :: start_time, end_time

      macrogrid_size_x = new_macrogrid_size_x
      macrogrid_size_y = new_macrogrid_size_y
      subgrid_size = new_subgrid_size

      interface_size = (macrogrid_size_x*(macrogrid_size_y-1) + macrogrid_size_y*(macrogrid_size_x-1))*(subgrid_size-2)
      dx = 1.0d0/dble(macrogrid_size_x*subgrid_size-(macrogrid_size_x-1)-1)
      dy = 1.0d0/dble(macrogrid_size_y*subgrid_size-(macrogrid_size_y-1)-1)

      macrogrid => new_macrogrid
      subgrid_solver_method => new_subgrid_solver_method

      if (use_openmp) then
         subgrids_computer => compute_subgrids_openmp
         start_time = omp_get_wtime()
         call new_macrogrid_solver_method()
         end_time = omp_get_wtime()
      else
         subgrids_computer => compute_subgrids
         call cpu_time(start_time)
         call new_macrogrid_solver_method()
         call cpu_time(end_time)
      end if

      time = end_time - start_time

   end subroutine run_macrogrid_solver
   ! Возвращает результаты последнего запуска макрорешателя.
   ! Аргументы:
   ! res_time: время выполнения в секундах.
   ! res_iter: число выполненных макроитераций.
   subroutine get_macrogrid_solver_results(res_time, res_iter)
      implicit none
      real*8, intent(out) :: res_time
      integer, intent(out) :: res_iter
      res_time = time
      res_iter = iter
   end subroutine get_macrogrid_solver_results
   ! Последовательно запускает решатель на каждой подсетке макросетки.
   subroutine compute_subgrids()
      implicit none
      integer :: iX, iY

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver_method(macrogrid(iX,iY,:,:), subgrid_size)
         end do
      end do

   end subroutine compute_subgrids
   ! Параллельно запускает решатель на каждой подсетке макросетки через OpenMP.
   subroutine compute_subgrids_openmp()
      implicit none
      integer :: iX, iY

      !$OMP PARALLEL DO PRIVATE(iX, iY) COLLAPSE(2) SCHEDULE(static) &
      !$OMP DEFAULT(NONE) SHARED(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver_method(macrogrid(iX,iY,:,:), subgrid_size)
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine compute_subgrids_openmp
   ! Выполняет простой итерационный метод согласования интерфейсов между подсетками.
   subroutine simple_iteration()
      implicit none
      real*8 :: old_vec(interface_size), new_vec(interface_size)
      real*8 :: interface_error
      integer :: iX, iY, k
      real*8 :: new_val, old_value

      do iter = 1, max_iter

         interface_error = 0.0d0
         call subgrids_computer()

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (4.0d0*(macrogrid(iX, iY, k, subgrid_size-1) + macrogrid(iX, iY+1, k, 2)) - &
                     macrogrid(iX, iY, k, subgrid_size-2) - macrogrid(iX, iY+1, k, 3))/6.0d0

                  macrogrid(iX, iY, k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (4.0d0*(macrogrid(iX, iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) - &
                     macrogrid(iX,iY, subgrid_size-2, k) - macrogrid(iX+1, iY, 3, k))/6.d0

                  macrogrid(iX, iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1, k) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         if (interface_error < eps) then
            exit
         end if

      end do

      call compute_intersection_nodes()

   end subroutine simple_iteration
   ! Выполняет простой итерационный метод при одном шаге решателя подсетки на макроитерацию.
   subroutine simple_iteration_one_iter()
      implicit none
      integer :: max_subgrid_iter

      call get_subgrid_solver_settings(new_max_iter = max_subgrid_iter)
      call set_subgrid_solver_settings(new_max_iter = 1)
      call simple_iteration()
      call set_subgrid_solver_settings(new_max_iter = max_subgrid_iter)

   end subroutine simple_iteration_one_iter
   ! Выполняет метод верхней релаксации с автоматической подстройкой параметра omega.
   subroutine sor()
      implicit none
      real*8 :: old_vec(interface_size), new_vec(interface_size)
      real*8 :: interface_error
      integer :: iX, iY, k
      real*8 :: new_val, old_value, norm1, norm2, pz

      omega = 1.0d0

      do iter = 1, max_iter

         interface_error = 0.0d0
         norm2 = 0.0d0
         call subgrids_computer()

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (4.0d0*(macrogrid(iX, iY, k, subgrid_size-1) + macrogrid(iX, iY+1, k, 2)) - &
                     macrogrid(iX, iY, k, subgrid_size-2) - macrogrid(iX, iY+1, k, 3))/6.0d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1) = new_val

                  norm2 = norm2 + dabs(new_val - old_value)
                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (4.0d0*(macrogrid(iX, iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) - &
                     macrogrid(iX,iY, subgrid_size-2, k) - macrogrid(iX+1, iY, 3, k))/6.d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1, k) = new_val

                  norm2 = norm2 + dabs(new_val - old_value)
                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         if (iter > 1) then
            pz = dsqrt(norm2)/dsqrt(norm1)
            if (pz < 1.0d0) then
               omega = 2.0d0 / (1.0d0 + dsqrt(1.0d0 - pz))
            end if
         end if

         norm1 = norm2

         if (interface_error < eps) then
            exit
         end if

      end do

      call compute_intersection_nodes()

   end subroutine sor
   ! Выполняет метод верхней релаксации с фиксированным значением omega.
   subroutine sor_fixed_omega()
      implicit none
      real*8 :: old_vec(interface_size), new_vec(interface_size)
      real*8 :: interface_error
      integer :: iX, iY, k
      real*8 :: new_val, old_value

      do iter = 1, max_iter

         interface_error = 0.0d0
         call subgrids_computer()

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (4.0d0*(macrogrid(iX, iY, k, subgrid_size-1) + macrogrid(iX, iY+1, k, 2)) - &
                     macrogrid(iX, iY, k, subgrid_size-2) - macrogrid(iX, iY+1, k, 3))/6.0d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (4.0d0*(macrogrid(iX, iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) - &
                     macrogrid(iX,iY, subgrid_size-2, k) - macrogrid(iX+1, iY, 3, k))/6.d0

                  new_val = new_val * omega + (1.0d0-omega) * old_value

                  macrogrid(iX, iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1, k) = new_val

                  interface_error = interface_error + dabs(new_val - old_value)

               end do
            end do
         end do

         if (interface_error < eps) then
            exit
         end if

      end do

      call compute_intersection_nodes()

   end subroutine sor_fixed_omega
   ! Выполняет SOR с фиксированным omega и одним шагом решателя подсетки на макроитерацию.
   subroutine sor_fixed_omega_one_iter()
      implicit none
      integer :: max_subgrid_iter

      call get_subgrid_solver_settings(new_max_iter = max_subgrid_iter)
      call set_subgrid_solver_settings(new_max_iter = 1)
      call sor_fixed_omega()
      call set_subgrid_solver_settings(new_max_iter = max_subgrid_iter)

   end subroutine sor_fixed_omega_one_iter
   ! Решает интерфейсную задачу методом сопряженных невязок.
   subroutine conjugate_residuals()
      implicit none
      real(8), dimension(interface_size) :: b, u_new, u_old, r_old, r_new, p_old, p_new
      real(8), dimension(interface_size) :: Ar_old, Ar_new, Ap_old, Ap_new
      real(8) :: alpha, beta, temp0, temp1, old_value, new_val, interface_error
      integer :: i

      call subgrids_computer()
      call s(b)

      u_old = 0.0d0

      call Cv(u_old, b, Ar_old)
      r_new = - b - Ar_old

      p_new = r_new

      call Cv(r_new, b, Ar_new)
      temp0 = dot_product(Ar_new, r_new)
      temp1 = dot_product(Ar_new, Ar_new)
      alpha = temp0 / temp1

      r_old = r_new - alpha * Ar_new

      call Cv(r_old, b, Ar_old)
      temp0 = dot_product(Ar_old, r_old)
      temp1 = dot_product(Ar_new, r_new)
      beta = temp0 / temp1

      p_old = r_old + beta * p_new
      Ap_old = Ar_old + beta * Ar_new
      u_old = alpha * p_new

      do iter = 1, max_iter

         temp0 = dot_product(Ar_old, r_old)
         temp1 = dot_product(Ap_old, Ap_old)
         alpha = temp0 / temp1

         u_new = u_old + alpha * p_old

         interface_error = 0.0d0
         do i = 1, interface_size
            interface_error = interface_error + dabs(u_new(i) - u_old(i))
         end do
         if (interface_error < eps) then
            exit
         end if

         r_new = r_old - alpha * Ap_old

         call Cv(r_new, b, Ar_new)
         temp0 = dot_product(Ar_new, r_new)
         temp1 = dot_product(Ar_old, r_old)
         beta = temp0 / temp1

         p_old = r_new + beta * p_old
         Ap_old = Ar_new + beta * Ap_old

         r_old = r_new
         Ar_old = Ar_new
         u_old = u_new

      end do

      call set_interface(u_new)
      call subgrids_computer()
      call compute_intersection_nodes()

   end subroutine conjugate_residuals
   ! Усредняет значения в узлах пересечения нескольких интерфейсов подсеток.
   subroutine compute_intersection_nodes()
      implicit none
      integer :: iX, iY
      real*8 :: old_value, new_val

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y - 1

            old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)

            new_val = (macrogrid(iX, iY, subgrid_size-1, subgrid_size) + &
               macrogrid(iX+1, iY, 2, subgrid_size) + &
               macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
               macrogrid(iX, iY+1, subgrid_size, 2 ))/4.d0

            macrogrid(iX, iY, subgrid_size, subgrid_size) = new_val
            macrogrid(iX+1, iY, 1, subgrid_size) = new_val
            macrogrid(iX, iY+1, subgrid_size, 1) = new_val
            macrogrid(iX+1, iY+1, 1, 1) = new_val

         end do
      end do

   end subroutine compute_intersection_nodes
   ! Записывает вектор интерфейсных значений обратно в макросетку.
   ! Аргументы:
   ! vec: линейный вектор значений на внутренних интерфейсах подсеток.
   subroutine set_interface(vec)
      implicit none
      real*8, intent(in) :: vec(interface_size)
      integer :: iX, iY, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               macrogrid(iX, iY, k, subgrid_size) = vec(i)
               macrogrid(iX, iY+1, k, 1) = vec(i)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               macrogrid(iX, iY, subgrid_size, k) = vec(i)
               macrogrid(iX+1, iY, 1, k) = vec(i)
               i = i + 1

            end do
         end do
      end do

   end subroutine set_interface
   ! Собирает текущие интерфейсные значения макросетки в линейный вектор.
   ! Аргументы:
   ! s_vec: выходной вектор интерфейсных значений.
   subroutine s(s_vec)
      implicit none
      real*8,  intent(out) :: s_vec(interface_size)

      integer :: iX, iY, i1, j1, k, i

      i = 1
      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               s_vec(i) = (3.0d0*(macrogrid(iX, iY, k, subgrid_size) + macrogrid(iX, iY+1, k, 1)) - &
                  4.0d0*(macrogrid(iX, iY, k, subgrid_size-1 ) + macrogrid(iX, iY+1, k, 2 )) + &
                  macrogrid(iX, iY, k, subgrid_size-2) + macrogrid(iX, iY+1, k, 3 ))/(2.0d0*dy)
               i = i + 1

            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               s_vec(i) = (3.0d0*(macrogrid(iX, iY, subgrid_size, k) + macrogrid(iX+1, iY, 1, k)) - &
                  4.0d0*(macrogrid(iX,iY, subgrid_size-1, k) + macrogrid(iX+1, iY, 2, k)) + &
                  macrogrid(iX,iY, subgrid_size-2, k) + macrogrid(iX+1, iY, 3, k))/(2.0d0*dx)
               i = i + 1

            end do
         end do
      end do

   end subroutine s
   ! Вычисляет действие интерфейсного оператора на вектор с учетом правой части.
   ! Аргументы:
   ! v: входной вектор интерфейсных поправок.
   ! b: базовый вектор интерфейсных значений.
   ! Cv_vec: результат применения оператора к вектору v.
   subroutine Cv(v, b, Cv_vec)
      implicit none
      real*8, intent(in) :: v(interface_size), b(interface_size)
      real*8, intent(out) :: Cv_vec(interface_size)

      real*8 :: subgrid_error

      call set_interface(v)
      call subgrids_computer()
      call s(Cv_vec)

      Cv_vec = Cv_vec - b

   end subroutine Cv

end module macrogrid_solver
