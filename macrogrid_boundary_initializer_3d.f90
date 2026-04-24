module macrogrid_boundary_initializer_3d
   implicit none
   private

   public :: initialize_macrogrid_constant_boundary_3d
   public :: compute_macrogrid_constant_boundary_error_3d

contains

   ! Инициализирует границу трехмерной макросетки значением 1.0.
   ! Аргументы:
   ! macrogrid: шестимерный массив значений по блокам и локальным узлам внутри блока.
   ! nx: число блоков по оси X.
   ! ny: число блоков по оси Y.
   ! nz: число блоков по оси Z.
   ! s: размер одной кубической подсетки по одной координате.
   subroutine initialize_macrogrid_constant_boundary_3d(macrogrid, nx, ny, nz, s)
      implicit none
      integer, intent(in) :: nx, ny, nz, s
      real*8,  intent(inout) :: macrogrid(nx, ny, nz, s, s, s)

      integer :: g_nx, g_ny, g_nz
      integer :: ix, iy, iz, i, j, k, lx, ly, lz

      g_nx = nx * s - (nx - 1)
      g_ny = ny * s - (ny - 1)
      g_nz = nz * s - (nz - 1)

      do ix = 1, nx
         do iy = 1, ny
            do iz = 1, nz
               do i = 1, s
                  do j = 1, s
                     do k = 1, s
                        lx = i + (ix-1)*s - (ix-1)
                        ly = j + (iy-1)*s - (iy-1)
                        lz = k + (iz-1)*s - (iz-1)

                        if (lx == 1 .or. lx == g_nx .or. &
                            ly == 1 .or. ly == g_ny .or. &
                            lz == 1 .or. lz == g_nz) then
                           macrogrid(ix, iy, iz, i, j, k) = 1.0d0
                        else
                           macrogrid(ix, iy, iz, i, j, k) = 0.0d0
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
      
   end subroutine initialize_macrogrid_constant_boundary_3d

   ! Вычисляет максимальное отклонение значений макросетки от 1.0.
   ! Аргументы:
   ! macrogrid: шестимерный массив текущего численного решения.
   ! nx: число блоков по оси X.
   ! ny: число блоков по оси Y.
   ! nz: число блоков по оси Z.
   ! s: размер одной кубической подсетки по одной координате.
   ! error: возвращаемая максимальная абсолютная ошибка.
   subroutine compute_macrogrid_constant_boundary_error_3d(macrogrid, nx, ny, nz, s, error)
      implicit none
      integer, intent(in) :: nx, ny, nz, s
      real*8,  intent(in) :: macrogrid(nx, ny, nz, s, s, s)
      real*8,  intent(out) :: error

      integer :: ix, iy, iz, i, j, k
      real*8  :: diff

      error = 0.0d0

      do ix = 1, nx
         do iy = 1, ny
            do iz = 1, nz
               do i = 1, s
                  do j = 1, s
                     do k = 1, s
                        diff = abs(macrogrid(ix,iy,iz,i,j,k) - 1.0d0)
                        if (diff > error) error = diff
                     end do
                  end do
               end do
            end do
         end do
      end do

   end subroutine compute_macrogrid_constant_boundary_error_3d

end module macrogrid_boundary_initializer_3d
