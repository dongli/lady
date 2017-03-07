module quad_points

  use commons
  use shape_function

  implicit none

  private

  public quad_points_init
  public quad_points_final

  real, public, allocatable :: Wq(:)         ! Quadrature weights
  real, public, allocatable :: Fq(:)         ! Shape function values
  real, public, allocatable :: Yq(:,:)       ! Quadrature body coordinates
  real, public, allocatable :: Xq(:,:,:)     ! Quadrature space coordinates
  real, public, allocatable :: Tq(:,:)       ! Temperature on quadrature points
  real, public, allocatable :: Vq(:,:,:)     ! Velocity on quadrature points
  real, public, allocatable :: Rhoq(:,:,:,:) ! Mass density on quadrature points

contains

  subroutine quad_points_init()

    integer qi, di, offset, ai(num_dim)

    num_quad_point = num_quad_point_span ** num_dim

    allocate(Wq(num_quad_point))
    allocate(Fq(num_quad_point))
    allocate(Yq(num_dim,num_quad_point))
    allocate(Xq(num_dim,num_quad_point,num_parcel))
    allocate(Tq(num_quad_point,num_parcel))
    allocate(Vq(num_dim,num_quad_point,num_parcel))
    allocate(Rhoq(num_dim,num_dim,num_quad_point,num_parcel))

    Wq(:) = 1.0d0
    do qi = 1, num_quad_point
      offset = 0.0d0
      do di = num_dim, 1, -1
        if (di /= num_dim) then
          offset = offset + (ai(di+1) - 1) * num_quad_point_span ** di
        end if
        ai(di) = ceiling(real(qi - offset) / num_quad_point_span ** (di - 1))
        Yq(di,qi) = quad_point_nodes(ai(di))
        Wq(qi) = Wq(qi) * quad_point_weights(ai(di))
      end do
      call shape_function_eval(Yq(:,qi), Fq(qi))
    end do

  end subroutine quad_points_init

  subroutine quad_points_final()

    deallocate(Wq)
    deallocate(Yq)
    deallocate(Xq)
    deallocate(Tq)
    deallocate(Vq)
    deallocate(Rhoq)

  end subroutine quad_points_final

end module quad_points
