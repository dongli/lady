module commons

  implicit none

  integer :: num_dim = -1
  integer :: num_parcel = -1
  integer :: num_quad_point = -1 ! Quadrature point number per parcel

  character(30) :: domain_type = 'cartesian'

end module commons
