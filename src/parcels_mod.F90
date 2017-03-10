module parcels_mod

  use commons_mod

  implicit none

  private

  public parcels_init
  public parcels_final
  public parcels_get_space_coord
  public parcels_get_body_coord

  public parcels

  type parcels_type
    real, allocatable :: X(:,:)   ! Centroid coordinate
    real, allocatable :: H(:,:,:) ! Deformation matrix
    real, allocatable :: I(:,:,:) ! Inverse of deformation matrix
    real, allocatable :: M(:)     ! Mass
    real, allocatable :: U(:)     ! Internal energy
  end type parcels_type

  type(parcels_type) parcels

contains

  subroutine parcels_init()

    allocate(parcels%X(commons%num_dim,commons%num_parcel))
    allocate(parcels%H(commons%num_dim,commons%num_dim,commons%num_parcel))
    allocate(parcels%I(commons%num_dim,commons%num_dim,commons%num_parcel))
    allocate(parcels%M(commons%num_parcel))
    if (commons%run_mode /= 'advection') then
      allocate(parcels%U(commons%num_parcel))
    end if

  end subroutine parcels_init

  subroutine parcels_reorder()

  end subroutine parcels_reorder

  subroutine parcels_final()

    if (allocated(parcels%X)) deallocate(parcels%X)
    if (allocated(parcels%H)) deallocate(parcels%H)
    if (allocated(parcels%I)) deallocate(parcels%I)
    if (allocated(parcels%M)) deallocate(parcels%M)
    if (allocated(parcels%U)) deallocate(parcels%U)

  end subroutine parcels_final

  subroutine parcels_get_space_coord(pi, y, x)

    integer, intent(in) :: pi
    real, intent(in) :: y(:)
    real, intent(out) :: x(:)

    x(:) = parcels%X(:,pi) + matmul(parcels%H(:,:,pi), y(:))

  end subroutine parcels_get_space_coord

  subroutine parcels_get_body_coord(pi, x, y)

    integer, intent(in) :: pi
    real, intent(in) :: x(:)
    real, intent(out) :: y(:)

    y(:) = matmul(parcels%I(:,:,pi), x(:) - parcels%X(:,pi))

  end subroutine parcels_get_body_coord

end module parcels_mod
