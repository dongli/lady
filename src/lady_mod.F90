module lady_mod

  use commons_mod
  use parcels_mod
  use quad_points_mod

  use cartesian_domain_mod

  implicit none

  private

  public lady_init
  public lady_final

contains

  subroutine lady_init(run_mode, num_parcel)

    character(*), intent(in) :: run_mode
    integer, intent(in) :: num_parcel

    commons%run_mode = run_mode
    commons%num_parcel = num_parcel

    call parcels_init()
    call quad_points_init()

  end subroutine lady_init

  subroutine lady_final()

    call parcels_final()
    call quad_points_final()

  end subroutine lady_final

  subroutine lady_calc_quad_point_space_coords()

    integer pi, qi

    do pi = 1, commons%num_parcel
      do qi = 1, commons%num_quad_point
        call parcels_get_space_coord(pi, quad_points%Y(:,qi), quad_points%X(:,qi,pi))
      end do
    end do

  end subroutine lady_calc_quad_point_space_coords

end module lady_mod
