module shape_function_mod

  use commons_mod

  implicit none

  private

  public shape_function_init
  public shape_function_eval
  public shape_function_diff

  public shape_function

  type shape_function_type
    integer num_quad_point_span ! Quadrature point number along each direction.
    real J
    real C
    real, allocatable :: quad_point_nodes(:)
    real, allocatable :: quad_point_weights(:)
  end type shape_function_type

  type(shape_function_type) shape_function

contains

  subroutine shape_function_init()

    shape_function%num_quad_point_span = 5
    shape_function%J = 1.0 / 12.0
    shape_function%C = 4.0 / 3.0

    allocate(shape_function%quad_point_nodes(shape_function%num_quad_point_span))
    allocate(shape_function%quad_point_weights(shape_function%num_quad_point_span))

    shape_function%quad_point_nodes(1)   =  -2.0d0/3.0d0
    shape_function%quad_point_nodes(2)   =  -1.0d0/3.0d0
    shape_function%quad_point_nodes(3)   =         0.0d0
    shape_function%quad_point_nodes(4)   =   1.0d0/3.0d0
    shape_function%quad_point_nodes(5)   =   2.0d0/3.0d0
    shape_function%quad_point_weights(1) =  41.0d0/1280.0d0
    shape_function%quad_point_weights(2) = 316.0d0/1280.0d0
    shape_function%quad_point_weights(3) = 566.0d0/1280.0d0
    shape_function%quad_point_weights(4) = 316.0d0/1280.0d0
    shape_function%quad_point_weights(5) =  41.0d0/1280.0d0

  end subroutine shape_function_init

  subroutine shape_function_eval(y, f)

    real, intent(in) :: y(:) ! Body coordinate
    real, intent(out) :: f   ! Shape function value

    integer i

    f = 1.0d0
    do i = 1, commons%num_dim
      if (-1.0d0 <= y(i) .and. y(i) <= -0.5d0) then
        f = f * (2.0d0 * (1.0d0 + y(i))**3.0d0)
      else if (-0.5d0 <= y(i) .and. y(i) <= 0.0d0) then
        f = f * (1.0d0 - 6.0d0 * y(i)**2.0d0 * (1.0d0 + y(i)))
      else if (0.0d0 <= y(i) .and. y(i) <= 0.5d0) then
        f = f * (1.0d0 - 6.0d0 * y(i)**2.0d0 * (1.0d0 - y(i)))
      else if (0.5d0 <= y(i) .and. y(i) <= 1.0d0) then
        f = f * (2.0d0 * (1.0d0 - y(i))**3.0d0)
      end if
    end do
    f = f * shape_function%C

  end subroutine shape_function_eval

  subroutine shape_function_diff(y, df)

    real, intent(in) :: y(:)   ! Body coordinate
    real, intent(out) :: df(:) ! Shape function derivatives

    integer i, j

    df(:) = 1.0d0
    do j = 1, commons%num_dim
      do i = 1, commons%num_dim
        if (df(j) == 0.0d0) cycle
        if (i == j) then
          if (-1.0d0 <= y(i) .and. y(i) <= -0.5d0) then
            df(j) = df(j) *    6.0d0 * (1.0d0 + y(i))**2.0d0
          else if (-0.5d0 <= y(i) .and. y(i) <= 0.0d0) then
            df(j) = df(j) * (- 12.0d0 * y(i) - 18.0d0 * y(i)**2.0d0)
          else if (0.0d0 <= y(i) .and. y(i) <= 0.5d0) then
            df(j) = df(j) * (- 12.0d0 * y(i) + 18.0d0 * y(i)**2.0d0)
          else if (0.5d0 <= y(i) .and. y(i) <= 1.0d0) then
            df(j) = df(j) * (- 6.0d0 * (1.0d0 - y(i))**2.0d0)
          else
            df(j) = 0.0d0
          end if
        else
          if (-1.0d0 <= y(i) .and. y(i) <= -0.5d0) then
            df(j) = df(j) * (2.0d0 * (1.0d0 + y(i))**3.0d0)
          else if (-0.5d0 <= y(i) .and. y(i) <= 0.0d0) then
            df(j) = df(j) * (1.0d0 - 6.0d0 * y(i)**2.0d0 * (1.0d0 + y(i)))
          else if (0.0d0 <= y(i) .and. y(i) <= 0.5d0) then
            df(j) = df(j) * (1.0d0 - 6.0d0 * y(i)**2.0d0 * (1.0d0 - y(i)))
          else if (0.5d0 <= y(i) .and. y(i) <= 1.0d0) then
            df(j) = df(j) * (2.0d0 * (1.0d0 - y(i))**3.0d0)
          else
            df(j) = 0.0d0
          end if
        end if
      end do
    end do
    df(:) = df(:) * shape_function%C

  end subroutine shape_function_diff

end module shape_function_mod
