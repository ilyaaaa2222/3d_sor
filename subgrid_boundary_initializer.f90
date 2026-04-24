module subgrid_boundary_initializer
   implicit none
   private

   public :: initialize_subgrid_boundary
   public :: compute_subgrid_boundary_error

   ! Параметры аналитического решения для двумерной краевой задачи.
   ! R1, R2: внутренний и внешний радиусы эталонного логарифмического профиля.
   ! x_min, y_min: координаты левого нижнего угла расчетной области.
   ! len: длина стороны квадратной области.
   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

contains

   ! Инициализирует границы двумерной подсетки аналитическим профилем.
   ! Аргументы:
   ! u: одномерное представление квадратной подсетки размера u_size x u_size.
   ! u_size: число узлов подсетки по одной координате.
   subroutine initialize_subgrid_boundary(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout) :: u(u_size*u_size)

      integer :: l0, l1, i

      do l1 = 1, u_size
         do l0 = 1, u_size

            i = l0 + u_size*(l1 - 1)

            if (l0.eq.1 .or. l0.eq.u_size .or. l1.eq.1 .or. l1.eq.u_size) then
               u(i) = log(sqrt((x_min + len*(l1 - 1)/u_size)**2 + &
                  (y_min + len*(l0 - 1)/u_size)**2)*R2/(R1*R1))/log(R2/R1)
            else
               u(i) = 0.0d0
            end if

         end do
      end do

   end subroutine initialize_subgrid_boundary

   ! Вычисляет максимальную ошибку подсетки относительно аналитического профиля.
   ! Аргументы:
   ! u: одномерное представление квадратной подсетки размера u_size x u_size.
   ! u_size: число узлов подсетки по одной координате.
   ! error: возвращаемая максимальная абсолютная ошибка.
   subroutine compute_subgrid_boundary_error(u, u_size, error)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout) :: u(u_size*u_size)
      real*8,  intent(out) :: error

      integer :: l0, l1, i

      error = 0.0d0

      do l1 = 1, u_size
         do l0 = 1, u_size

            i = l0 + u_size*(l1 - 1)
            error = max(error, abs(u(i) - log(sqrt((x_min + len*(l1 - 1)/u_size)**2 + &
               (y_min + len*(l0 - 1)/u_size)**2)*R2/(R1*R1))/log(R2/R1)))

         end do
      end do

   end subroutine compute_subgrid_boundary_error

end module subgrid_boundary_initializer
