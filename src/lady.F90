module lady

  use parcels
  use quad_points

contains

  subroutine lady_init()

    call parcels_init()
    call quad_points_init()

  end subroutine lady_init

end module lady
