module commons_mod

  implicit none

  type commons_type
    integer :: num_dim = -1
    integer :: num_parcel = -1
    integer :: num_quad_point = -1 ! Quadrature point number per parcel

    character(30) :: domain_type = 'cartesian'
    character(30) :: run_mode = '' ! Or 'shallow_water', 'barotropic', 'baroclinic', or 'advection'
  end type commons_type

  type(commons_type) commons

end module commons_mod
