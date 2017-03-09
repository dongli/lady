module lady

  use commons
  use parcels
  use quad_points

  use cartesian_domain

  implicit none

  private

  public lady_init
  public lady_final

  character(30) domain

contains

  subroutine lady_init(run_mode_, num_parcel_)

    character(*), intent(in) :: run_mode_
    integer, intent(in) :: num_parcel_

    run_mode = run_mode_
    num_parcel = num_parcel_

    call parcels_init()
    call quad_points_init()

  end subroutine lady_init

  subroutine lady_final()

    call parcels_final()
    call quad_points_final()

  end subroutine lady_final

  subroutine lady_calc_quad_point_space_coords()

    integer pi, qi

    do pi = 1, num_parcel
      do qi = 1, num_quad_point
        call parcels_get_space_coord(pi, Yq(:,qi), Xq(:,qi,pi))
      end do
    end do

  end subroutine lady_calc_quad_point_space_coords

end module lady
