module cartesian_domain

  use commons, only: domain_type, num_dim

  implicit none

  private

  public cartesian_domain_init
  public cartesian_domain_final

  integer, public, parameter :: BT_OPEN = 1
  integer, public, parameter :: BT_WALL = 2
  integer, public, parameter :: BT_PERIODIC = 3

  real :: axis_x_range(2) = [0, 1]
  integer :: axis_x_btype(2) = [BT_OPEN, BT_OPEN]
  real :: axis_y_range(2) = [0, 1]
  integer :: axis_y_btype(2) = [BT_OPEN, BT_OPEN]
  real :: axis_z_range(2) = [0, 1]
  integer :: axis_z_btype(2) = [BT_WALL, BT_OPEN]

contains

  subroutine cartesian_domain_init(num_dim_, &
      axis_x_range_, axis_x_btype_, &
      axis_y_range_, axis_y_btype_, &
      axis_z_range_, axis_z_btype_)

    integer, intent(in) :: num_dim_
    real, intent(in), optional :: axis_x_range_(2)
    integer, intent(in), optional :: axis_x_btype_(2)
    real, intent(in), optional :: axis_y_range_(2)
    integer, intent(in), optional :: axis_y_btype_(2)
    real, intent(in), optional :: axis_z_range_(2)
    integer, intent(in), optional :: axis_z_btype_(2)

    domain_type = 'cartesian'
    num_dim = num_dim_

    if (present(axis_x_range_)) axis_x_range = axis_x_range_
    if (present(axis_x_btype_)) axis_x_btype = axis_x_btype_
    if (present(axis_y_range_)) axis_y_range = axis_y_range_
    if (present(axis_y_btype_)) axis_y_btype = axis_y_btype_
    if (present(axis_z_range_)) axis_z_range = axis_z_range_
    if (present(axis_z_btype_)) axis_z_btype = axis_z_btype_

  end subroutine cartesian_domain_init

  subroutine cartesian_domain_final()

  end subroutine cartesian_domain_final

end module cartesian_domain
