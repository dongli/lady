module quad_points_mod

  use commons_mod
  use shape_function_mod

  implicit none

  private

  public quad_points_init
  public quad_points_final

  public quad_points

  type quad_points_type
    real, allocatable :: W(:)     ! Quadrature weights
    real, allocatable :: F(:)     ! Shape function values
    real, allocatable :: Y(:,:)   ! Quadrature body coordinates
    real, allocatable :: X(:,:,:) ! Quadrature space coordinates
    real, allocatable :: T(:,:)   ! Temperature on quadrature points
    real, allocatable :: V(:,:,:) ! Velocity on quadrature points
    real, allocatable :: rho(:,:) ! Mass density on quadrature points
  end type quad_points_type

  type(quad_points_type) quad_points

contains

  subroutine quad_points_init()

    integer qi, di, offset, ai(commons%num_dim)

    if (commons%run_mode == 'advection') return

    commons%num_quad_point = shape_function%num_quad_point_span ** commons%num_dim

    allocate(quad_points%W(commons%num_quad_point))
    allocate(quad_points%F(commons%num_quad_point))
    allocate(quad_points%Y(commons%num_dim,commons%num_quad_point))
    allocate(quad_points%X(commons%num_dim,commons%num_quad_point,commons%num_parcel))
    allocate(quad_points%T(commons%num_quad_point,commons%num_parcel))
    allocate(quad_points%V(commons%num_dim,commons%num_quad_point,commons%num_parcel))
    allocate(quad_points%rho(commons%num_quad_point,commons%num_parcel))

    quad_points%W(:) = 1.0d0
    do qi = 1, commons%num_quad_point
      offset = 0.0d0
      do di = commons%num_dim, 1, -1
        if (di /= commons%num_dim) then
          offset = offset + (ai(di+1) - 1) * shape_function%num_quad_point_span ** di
        end if
        ai(di) = ceiling(real(qi - offset) / shape_function%num_quad_point_span ** (di - 1))
        quad_points%Y(di,qi) = shape_function%quad_point_nodes(ai(di))
        quad_points%W(qi) = quad_points%W(qi) * shape_function%quad_point_weights(ai(di))
      end do
      call shape_function_eval(quad_points%Y(:,qi), quad_points%F(qi))
    end do

  end subroutine quad_points_init

  subroutine quad_points_final()

    if (allocated(quad_points%W)) deallocate(quad_points%W)
    if (allocated(quad_points%Y)) deallocate(quad_points%Y)
    if (allocated(quad_points%X)) deallocate(quad_points%X)
    if (allocated(quad_points%T)) deallocate(quad_points%T)
    if (allocated(quad_points%V)) deallocate(quad_points%V)
    if (allocated(quad_points%rho)) deallocate(quad_points%rho)

  end subroutine quad_points_final

end module quad_points_mod
