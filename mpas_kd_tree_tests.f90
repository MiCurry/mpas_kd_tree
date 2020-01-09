program mpas_kd_tester

   use mpas_kd_tree, only : kdnode 
   use mpas_kd_tree, only : mpas_kd_construct, mpas_kd_search
   use mpas_kd_tree, only : mpas_kd_free

   implicit none

   integer :: argc
   character (len=255) :: argv
   character :: c

   integer :: npoints
   integer :: ndims

   integer :: lrange
   integer :: urange

   integer :: ntests
   integer :: npoints_l
   integer :: npoints_u
   integer :: ndims_l
   integer :: ndims_u

   namelist /testVars/ npoints, ndims, lrange, urange, ntests, npoints_l, npoints_u, ndims_l, ndims_u

   open(42,file='namelist.input', status='old',form='formatted',action='read')
   read(42, testVars)
   close(43)

   call test1(npoints, ndims, lrange, urange)

contains

! Test create, search and free 1, 2, 3 and 7 elements
subroutine test1(npoints, ndims, lrange, urange)

   implicit none

   integer, intent(in), value :: npoints
   integer, intent(in), value :: ndims
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange

   real, dimension(:,:), pointer :: points
   type(kdnode), dimension(:), pointer :: nodes
   type(kdnode), pointer :: tree => null()
   type(kdnode), pointer :: res => null()


   write(0,*) "MPAS_KD_TREE Tests Started"

   !
   ! Testing mpas_kd_construct with 1 node
   !
   write(0,*) "Testing mpas_kd_construct with 1 node"

   allocate(nodes(1))
   allocate(nodes(1) % point(3))
   nodes(1) % point(:) = (/1.0, 1.0, 1.0/)
   nodes(1) % cell = 1

   tree => mpas_kd_construct(nodes, 3)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
   endif

   ! Root node should contain our point
   if (all(tree % point(:) /= (/1.0, 1.0, 1.0/))) then
      write(0,*) "ERROR: Root node did not contain the correct point"
   endif

   call mpas_kd_search(tree, (/1.0, 1.0, 1.0/), res)
   if (all(res % point(:) /= (/1.0, 1.0, 1.0/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % cell /= 1) then
      write(0,*) "ERROR: Cell did not match the expected cell"
   endif

   call mpas_kd_search(tree, (/1000.0, 1000.0, 1000.0/), res)
   if (all(res % point(:) /= (/1.0, 1.0, 1.0/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % cell /= 1) then
      write(0,*) "ERROR: Cell did not match the expected cell"
   endif

   call mpas_kd_free(tree)
   nullify(tree)
   deallocate(nodes)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
   endif

   ! 
   ! Testing mpas_kd_construct with 2 nodes
   !
   write(0,*) "Testing mpas_kd_construct with 2 node"
   allocate(nodes(2))

   allocate(nodes(1) % point(3))
   nodes(1) % point(:) = (/1.0, 1.0, 1.0/)
   nodes(1) % cell = 1

   allocate(nodes(2) % point(3))
   nodes(2) % point(:) = (/-1.0, -1.0, -1.0/)
   nodes(2) % cell = 2

   tree => mpas_kd_construct(nodes, 3)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
   endif

   call mpas_kd_search(tree, (/1.0, 1.0, 1.0/), res)
   if (all(res % point(:) /= (/1.0, 1.0, 1.0/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % cell /= 1) then
      write(0,*) "ERROR: Cell did not match the expected cell"
   endif

   call mpas_kd_search(tree, (/-1.0, -1.0, -1.0/), res)
   if (all(res % point(:) /= (/-1.0, -1.0, -1.0/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % cell /= 2) then
      write(0,*) "ERROR: Cell did not match the expected cell"
   endif

   call mpas_kd_free(tree)
   deallocate(nodes)
   nullify(tree)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
   endif

   write(0,*) "End"

end subroutine test1

end program mpas_kd_tester
