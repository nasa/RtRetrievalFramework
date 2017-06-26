module radiant_keywords

  implicit none

  integer, parameter :: RADIANT_BRDF_KERNEL_NONE      = -1
  integer, parameter :: RADIANT_BRDF_KERNEL_LAMBERT   = 0
  integer, parameter :: RADIANT_BRDF_KERNEL_ROSSTHIN  = 1
  integer, parameter :: RADIANT_BRDF_KERNEL_ROSSTHICK = 2
  integer, parameter :: RADIANT_BRDF_KERNEL_LISPARSE  = 3
  integer, parameter :: RADIANT_BRDF_KERNEL_LIDENSE   = 4
  integer, parameter :: RADIANT_BRDF_KERNEL_HAPKE     = 5
  integer, parameter :: RADIANT_BRDF_KERNEL_ROUJEAN   = 6
  integer, parameter :: RADIANT_BRDF_KERNEL_RAHMAN    = 7
  integer, parameter :: RADIANT_BRDF_KERNEL_COXMUNK   = 8
  integer, parameter :: RADIANT_BRDF_KERNEL_RHERMAN   = 9
  integer, parameter :: RADIANT_BRDF_KERNEL_BREON     = 10

  integer, parameter :: RADIANT_QUADRATURE_GAUSS        = 1
  integer, parameter :: RADIANT_QUADRATURE_DOUBLE_GAUSS = 2
  integer, parameter :: RADIANT_QUADRATURE_LOBATTO      = 3
  
  integer, parameter :: RADIANT_TASK_ALLOCATE_ALL      = 1
  integer, parameter :: RADIANT_TASK_DEALLOCATE_ALL    = 2

  integer, parameter :: RADIANT_TASK_ALLOCATE_STREAM   = 3
  integer, parameter :: RADIANT_TASK_DEALLOCATE_STREAM = 4

end module radiant_keywords
