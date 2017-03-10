module cartesian_mesh_mod

  use commons_mod

  implicit none

  private

  public cartesian_mesh_init
  public cartesian_mesh_final

  integer, public, parameter :: FULL = 1
  integer, public, parameter :: HALF = 2

  type cartesian_mesh_type
    integer, public :: nx(2)
    integer, public :: ny(2)
    integer, public :: nz(2)

    real, public, allocatable :: x(:,:) ! Grid X coordinates
    real, public, allocatable :: y(:,:) ! Grid Y coordinates
    real, public, allocatable :: z(:,:) ! Grid Z coordinates
  end type cartesian_mesh_type

  type(cartesian_mesh_type) mesh

contains

  subroutine cartesian_mesh_init(nx, ny, nz)

    integer, optional, intent(in) :: nx(2)
    integer, optional, intent(in) :: ny(2)
    integer, optional, intent(in) :: nz(2)

    if (present(nx)) then
      mesh%nx = nx
      allocate(mesh%x(maxval(nx),2))
    end if
    if (present(ny)) then
      mesh%ny = ny
      allocate(mesh%y(maxval(ny),2))
    end if
    if (present(nz)) then
      mesh%nz = nz
      allocate(mesh%z(maxval(nz),2))
    end if

  end subroutine cartesian_mesh_init

  subroutine cartesian_mesh_final()

    if (allocated(mesh%x)) deallocate(mesh%x)
    if (allocated(mesh%y)) deallocate(mesh%y)
    if (allocated(mesh%z)) deallocate(mesh%z)

  end subroutine cartesian_mesh_final

end module cartesian_mesh_mod
