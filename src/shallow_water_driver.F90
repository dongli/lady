program shallow_water_driver

  use lady_mod
  use cartesian_domain_mod

  implicit none

  call cartesian_domain_init(num_dim=2, &
    axis_x_range=[0.0, 100.0], axis_x_btype=[BT_OPEN, BT_OPEN], &
    axis_y_range=[0.0, 100.0], axis_y_btype=[BT_OPEN, BT_OPEN])
  call lady_init(run_mode='shallow_water', num_parcel=10000)

  call cartesian_domain_final()
  call lady_final()

end program shallow_water_driver
