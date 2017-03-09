program shallow_water_driver

  use lady
  use cartesian_domain

  implicit none

  call cartesian_domain_init(num_dim_=2, &
    axis_x_range_=[0.0, 100.0], axis_x_btype_=[BT_OPEN, BT_OPEN], &
    axis_y_range_=[0.0, 100.0], axis_y_btype_=[BT_OPEN, BT_OPEN])
  call lady_init(run_mode_='shallow_water', num_parcel_=10000)

  call cartesian_domain_final()
  call lady_final()

end program shallow_water_driver
