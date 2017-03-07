module parcels

  use commons

  implicit none

  private

  public parcels_init
  public parcels_final
  public parcels_get_space_coord
  public parcels_get_body_coord

  real, allocatable :: Xp(:,:)   ! Centroid coordinate
  real, allocatable :: Hp(:,:,:) ! Deformation matrix
  real, allocatable :: Ip(:,:,:) ! Inverse of deformation matrix
  real, allocatable :: Mp(:)     ! Mass
  real, allocatable :: Up(:)     ! Internal energy

contains

  subroutine parcels_init()

    allocate(Xp(num_dim,num_parcel))
    allocate(Hp(num_dim,num_dim,num_parcel))
    allocate(Ip(num_dim,num_dim,num_parcel))
    allocate(Mp(num_parcel))
    allocate(Up(num_parcel))

  end subroutine parcels_init

  subroutine parcels_reorder()

  end subroutine parcels_reorder

  subroutine parcels_final()

    deallocate(Xp)
    deallocate(Hp)
    deallocate(Mp)
    deallocate(Up)

  end subroutine parcels_final

  subroutine parcels_get_space_coord(pi, y, x)

    integer, intent(in) :: pi
    real, intent(in) :: y(:)
    real, intent(out) :: x(:)

    x(:) = Xp(:,pi) + matmul(Hp(:,:,pi), y(:))

  end subroutine parcels_get_space_coord

  subroutine parcels_get_body_coord(pi, x, y)

    integer, intent(in) :: pi
    real, intent(in) :: x(:)
    real, intent(out) :: y(:)

    y(:) = matmul(Ip(:,:,pi), x(:) - Xp(:,pi))

  end subroutine parcels_get_body_coord

end module parcels
