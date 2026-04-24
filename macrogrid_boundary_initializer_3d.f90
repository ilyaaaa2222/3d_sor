module macrogrid_boundary_initializer_3d
   implicit none
   private

   public :: initialize_macrogrid_constant_boundary_3d
   public :: compute_macrogrid_constant_boundary_error_3d

contains

   ! Инициализация: границы = 1.0, внутри = 0.0
   subroutine initialize_macrogrid_constant_boundary_3d(macrogrid, nx, ny, nz, s)
      implicit none
      integer, intent(in) :: nx, ny, nz, s
      ! Шестимерный массив: (блок_X, блок_Y, блок_Z, узел_i, узел_j, узел_k)
      real*8,  intent(inout) :: macrogrid(nx, ny, nz, s, s, s)

      integer :: g_nx, g_ny, g_nz
      integer :: ix, iy, iz, i, j, k, lx, ly, lz

      ! Глобальные размеры сетки (с учетом наложения границ подсеток)
      g_nx = nx * s - (nx - 1)
      g_ny = ny * s - (ny - 1)
      g_nz = nz * s - (nz - 1)

      do ix = 1, nx
         do iy = 1, ny
            do iz = 1, nz
               do i = 1, s
                  do j = 1, s
                     do k = 1, s
                        ! Считаем глобальный индекс узла в кубе
                        lx = i + (ix-1)*s - (ix-1)
                        ly = j + (iy-1)*s - (iy-1)
                        lz = k + (iz-1)*s - (iz-1)

                        ! Если узел лежит на любой из 6 внешних граней глобального куба
                        if (lx == 1 .or. lx == g_nx .or. &
                            ly == 1 .or. ly == g_ny .or. &
                            lz == 1 .or. lz == g_nz) then
                            
                            macrogrid(ix, iy, iz, i, j, k) = 1.0d0
                        else
                            ! Внутри куба — начальное приближение (ноль)
                            macrogrid(ix, iy, iz, i, j, k) = 0.0d0
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
      
   end subroutine initialize_macrogrid_constant_boundary_3d

   ! Расчет ошибки: насколько значения в сетке отклоняются от 1.0
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