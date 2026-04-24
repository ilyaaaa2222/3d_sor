module macrogrid_boundary_initializer
   implicit none
   private

   public :: initialize_macrogrid_boundary
   public :: initialize_macrogrid_constant_boundary
   public :: compute_macrogrid_boundary_error
   public :: compute_macrogrid_constant_boundary_error

   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

contains

   subroutine initialize_macrogrid_boundary(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

      integer :: global_size_x, global_size_y
      integer :: iX, iY, i1, j1, lX, lY
      real*8  :: boundary_val

      global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
      global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            do i1 = 1, subgrid_size
               do j1 = 1, subgrid_size
                  lX = i1 + (iX-1)*subgrid_size - (iX-1)
                  lY = j1 + (iY-1)*subgrid_size - (iY-1)

                  if (lX == 1 .or. lX == global_size_x .or. &
                     lY == 1 .or. lY == global_size_y) then
                     boundary_val = log( sqrt( (x_min + len*(lX-1)/(global_size_x-1))**2 + &
                        (y_min + len*(lY-1)/(global_size_y-1))**2 )*R2/(R1*R1) ) / &
                        log(R2/R1)

                     macrogrid(iX, iY, i1, j1) = boundary_val
                  else
                     macrogrid(iX, iY, i1, j1) = 0.0d0
                  end if
               end do
            end do
         end do
      end do
      
   end subroutine initialize_macrogrid_boundary

   subroutine initialize_macrogrid_constant_boundary(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
      real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

      integer :: global_size_x, global_size_y
      integer :: iX, iY, i1, j1, lX, lY

      global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
      global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            do i1 = 1, subgrid_size
               do j1 = 1, subgrid_size
                  lX = i1 + (iX-1)*subgrid_size - (iX-1)
                  lY = j1 + (iY-1)*subgrid_size - (iY-1)

                  if (lX == 1 .or. lX == global_size_x .or. &
                     lY == 1 .or. lY == global_size_y) then
                     macrogrid(iX, iY, i1, j1) = 1.0d0
                  else
                     macrogrid(iX, iY, i1, j1) = 0.0d0
                  end if
               end do
            end do
         end do
      end do

   end subroutine initialize_macrogrid_constant_boundary

   subroutine compute_macrogrid_boundary_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, error)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
      real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(out) :: error

      integer :: global_size_x, global_size_y
      integer :: iX, iY, i1, j1, lX, lY
      real*8  :: val_exact, diff

      global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
      global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

      error = 0.0d0

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            do i1 = 1, subgrid_size
               do j1 = 1, subgrid_size
                  lX = i1 + (iX-1)*subgrid_size - (iX-1)
                  lY = j1 + (iY-1)*subgrid_size - (iY-1)

                  val_exact = log( sqrt( (x_min + len*(lX-1)/(global_size_x-1))**2 + &
                     (y_min + len*(lY-1)/(global_size_y-1))**2 )*R2/(R1*R1) ) / &
                     log(R2/R1)

                  diff = abs(macrogrid(iX,iY,i1,j1) - val_exact)
                  if (diff > error) error = diff
               end do
            end do
         end do
      end do

   end subroutine compute_macrogrid_boundary_error

   subroutine compute_macrogrid_constant_boundary_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, error)
      implicit none
      integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
      real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
      real*8,  intent(out) :: error

      integer :: iX, iY, i1, j1
      real*8  :: diff

      error = 0.0d0

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            do i1 = 1, subgrid_size
               do j1 = 1, subgrid_size
                  diff = abs(macrogrid(iX,iY,i1,j1) - 1.0d0)
                  if (diff > error) error = diff
               end do
            end do
         end do
      end do

   end subroutine compute_macrogrid_constant_boundary_error

end module macrogrid_boundary_initializer
