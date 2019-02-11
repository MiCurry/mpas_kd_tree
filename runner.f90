program mpas_kd_tester

   use mpas_kd_tree, only : kdnode
   use mpas_kd_tree, only : mpas_kd_insert, mpas_kd_construct, mpas_kd_search, mpas_kd_remove
   use getoptf

   implicit none

   integer :: argc
   character (len=255) :: argv
   character :: c

   integer :: npoints
   integer :: ndims

   integer :: lrange
   integer :: urange

   namelist /testVars/ npoints, ndims, lrange, urange

   open(42,file='namelist.input', status='old',form='formatted',action='read')
   read(42, testVars)
   close(43)

   argc = command_argument_count()
   call get_command(argv)

   do while ( getopt(argc, argv, c, "12345") )
      select case (c)
         case ('1')
            write(0,*) "Launching test 1" 
            call test1(npoints, ndims, lrange, urange)
         case ('2')
            write(0,*) "Launching test 2"
            call test2(npoints, ndims, lrange, urange)
         case ('3')
            write(0,*) "Launching test 3"
         case ('4')
            write(0,*) "Launching test 4"
         case ('5')
            write(0,*) "Launching test 5"
      endselect
   enddo

contains


subroutine test1(npoints, ndims, lrange, urange)

   implicit none

   integer, intent(in), value :: npoints, ndims
   integer, intent(in), value :: lrange, urange

   real, dimension(ndims, npoints) :: r_points, rps, rpsr
   integer, dimension(ndims, npoints) :: i_points, ips, ipsr
   integer :: i

   type (kdnode), pointer :: iTree, rTree


end subroutine test1

subroutine test2(npoints, ndims, lrange, urange)

   implicit none
   integer, intent(in), value :: npoints, ndims
   integer, intent(in), value :: lrange, urange

   type(kdnode), pointer :: rTree => null()
   real, dimension(ndims, npoints) :: arry1, arry2
   real, dimension(ndims) :: result
   real :: min_d

   integer :: i, r, k, n
   
   min_d = huge(min_d)

   write(0,*) "Real test: "
   call random_number(arry1)
   call random_number(arry2)

   arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange
   
   write(0,*) ""
   !write(0,*) "arry 1: ", arry1
   write(0,*) ""
   !write(0,*) "arry 2: ", arry2
   write(0,*) ""

   rTree => mpas_kd_construct(arry1)

   write(0,*) ""

   if (associated(rTree)) then
      write(0,*) "Real Tree Node: ", rTree % data

      if (associated(rTree % left)) then
         write(0,*) "Left tree was associated"
         write(0,*) rTree % left % data
      else
         write(0,*) "Left tree was NOT associated"
      endif

      if (associated(rTree % right)) then
         write(0,*) "Right tree was associated"
         write(0,*) rTree % right % data
      else
         write(0,*) "Right tree was NOT associated"
      endif
   else
      write(0,*) "rTree was not associated :'("
   endif


   write(0,*) ""
   write(0,*) "Searching Test!"
   write(0,*) ""
   write(0,*) "Searching all points currently in the tree"
   write(0,*) ""

   do i = 1, size(arry1, dim=2)
      min_d = huge(min_d)
      call mpas_kd_search(rTree, arry1(:, i), result, min_d)
      if (.NOT. (all(result(:) == arry1(:,i)) .AND. min_d == 0.0)) then
         write(0,*) "That point wasn't in the tree, but its cloest point was: ", result, "and the dist: ", sqrt(min_d)
      endif
   enddo

   write(0,*) ""
   write(0,*) "Searching for randomly created points"
   write(0,*) ""

   do i = 1, size(arry2, dim=2)
      min_d = huge(min_d)
      call mpas_kd_search(rTree, arry2(:, i), result, min_d)
      if (.NOT. (all(result(:) == arry2(:,i)) .AND. min_d == 0.0)) then
         write(0,*) "That point wasn't in the tree, but its cloest point was: ", result, "and the dist: ", sqrt(min_d)
      endif
   enddo

end subroutine test2




end program mpas_kd_tester
