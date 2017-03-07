module shape_function

  use commons

  implicit none

  private

  public num_quad_point_span
  public shape_function_J
  public quad_point_nodes
  public quad_point_weights
  public shape_function_init
  public shape_function_eval
  public shape_function_diff

  integer, parameter :: num_quad_point_span = 5 ! Quadrature point number along each direction.
  real, parameter :: shape_function_J = 1.0d0 / 12.0d0
  real, parameter :: shape_function_C = 4.0d0 / 3.0d0

  real quad_point_nodes(num_quad_point_span)
  real quad_point_weights(num_quad_point_span)

contains

  subroutine shape_function_init()

    quad_point_nodes(1)   =  -2.0d0/3.0d0
    quad_point_nodes(2)   =  -1.0d0/3.0d0
    quad_point_nodes(3)   =         0.0d0
    quad_point_nodes(4)   =   1.0d0/3.0d0
    quad_point_nodes(5)   =   2.0d0/3.0d0
    quad_point_weights(1) =  41.0d0/1280.0d0
    quad_point_weights(2) = 316.0d0/1280.0d0
    quad_point_weights(3) = 566.0d0/1280.0d0
    quad_point_weights(4) = 316.0d0/1280.0d0
    quad_point_weights(5) =  41.0d0/1280.0d0

  end subroutine shape_function_init

  subroutine shape_function_eval(y, f)

    real, intent(in) :: y(:) ! Body coordinate
    real, intent(out) :: f   ! Shape function value

    integer i

    f = 1.0d0
    do i = 1, num_dim
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
    f = f * shape_function_C

  end subroutine shape_function_eval

  subroutine shape_function_diff(y, df)

    real, intent(in) :: y(:)   ! Body coordinate
    real, intent(out) :: df(:) ! Shape function derivatives

    integer i, j

    df(:) = 1.0d0
    do j = 1, num_dim
      do i = 1, num_dim
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
    df(:) = df(:) * shape_function_C

  end subroutine shape_function_diff

end module shape_function
