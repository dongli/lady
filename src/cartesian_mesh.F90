module cartesian_mesh

  use commons

  implicit none

  private

  public cartesian_mesh_init
  public cartesian_mesh_final

  integer, public, parameter :: FULL = 1
  integer, public, parameter :: HALF = 2

  integer, public :: nXg(2)
  integer, public :: nYg(2)
  integer, public :: nZg(2)

  real, public, allocatable :: Xg(:,:) ! Grid X coordinates
  real, public, allocatable :: Yg(:,:) ! Grid Y coordinates
  real, public, allocatable :: Zg(:,:) ! Grid Z coordinates

contains

  subroutine cartesian_mesh_init(nXg_, nYg_, nZg_)

    integer, optional, intent(in) :: nXg_(2)
    integer, optional, intent(in) :: nYg_(2)
    integer, optional, intent(in) :: nZg_(2)

    if (present(nXg_)) then
      nXg = nXg_
      allocate(Xg(maxval(nXg),2))
    end if
    if (present(nYg_)) then
      nYg = nYg_
      allocate(Yg(maxval(nYg),2))
    end if
    if (present(nZg_)) then
      nZg = nZg_
      allocate(Zg(maxval(nZg),2))
    end if

  end subroutine cartesian_mesh_init

  subroutine cartesian_mesh_final()

    if (allocated(Xg)) deallocate(Xg)
    if (allocated(Yg)) deallocate(Yg)
    if (allocated(Zg)) deallocate(Zg)

  end subroutine cartesian_mesh_final

end module cartesian_mesh
