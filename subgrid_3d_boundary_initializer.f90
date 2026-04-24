module subgrid_3d_boundary_initializer
   implicit none
   private

   public :: initialize_subgrid_3d_boundary
   public :: compute_subgrid_3d_boundary_error

   ! Параметры эталонной области, зарезервированные для задач с аналитической границей.
   ! Сейчас в трехмерном инициализаторе используется только константная граница 1.0.
   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

contains

   ! Инициализирует граничные узлы кубической подсетки значением 1.0.
   ! Аргументы:
   ! u: одномерное представление кубической подсетки размера u_size x u_size x u_size.
   ! u_size: число узлов подсетки по одной координате.
   subroutine initialize_subgrid_3d_boundary(u, u_size)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout), target :: u(u_size*u_size*u_size)

      integer :: l0, l1, l2, i

      do l2 = 1, u_size
         do l1 = 1, u_size
            do l0 = 1, u_size

               i = (l2-1)*u_size*u_size + (l1-1)*u_size + l0

               if (l2.eq.1 .or. l2.eq.u_size) then
                  u(i) = 1.0d0
               else if (l1.eq.1 .or. l1.eq.u_size) then
                  u(i) = 1.0d0
               else if (l0.eq.1 .or. l0.eq.u_size) then
                  u(i) = 1.0d0
               else
                  u(i) = 0.0d0
               end if

            end do
         end do
      end do
   end subroutine initialize_subgrid_3d_boundary

   ! Вычисляет максимальное отклонение значений подсетки от 1.0.
   ! Аргументы:
   ! u: одномерное представление кубической подсетки размера u_size x u_size x u_size.
   ! u_size: число узлов подсетки по одной координате.
   ! error: возвращаемая максимальная абсолютная ошибка.
   subroutine compute_subgrid_3d_boundary_error(u, u_size, error)
      implicit none
      integer, intent(in) :: u_size
      real*8,  intent(inout), target :: u(u_size*u_size*u_size)
      real*8,  intent(out) :: error

      integer :: l0, l1, l2, i

      error = 0.0d0

      do l2 = 1, u_size
         do l1 = 1, u_size
            do l0 = 1, u_size

               i = (l2-1)*u_size*u_size + (l1-1)*u_size + l0
               error = max(error, abs(u(i) - 1.0d0))

            end do
         end do
      end do

   end subroutine compute_subgrid_3d_boundary_error

end module subgrid_3d_boundary_initializer
