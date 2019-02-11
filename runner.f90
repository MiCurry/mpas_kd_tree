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
            call test2(ndims, lrange, urange)
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


   ! Integer Test
   call random_number(r_points)
   call random_number(rpsr)
   i_points = int(r_points(:,:) * (urange + 1 - lrange)) + lrange
   ipsr = int(rpsr(:,:) * (urange + 1 - lrange)) + lrange
   ips = i_points
   write(0,*) ""
   write(0,*) "=================="
   write(0,*) "Integer Test:"
   write(0,*) " Integer Shape: ", shape(i_points), " nElms: ", size(i_points)
   !write(0,*) " Integer Points: ", i_points
   write(0,*) ""

   write(0,*) "Constructing KD Tree of Integer Type..."
   iTree => mpas_kd_construct(i_points)
   write(0,*) "Finished mpas_kd_constrcut of an integer tree"

   write(0,*) "Searching for points..."
   i = mpas_kd_search(iTree, i_points(:,1), 1)
   write(0,*) "Finished mpas_kd_search of an integer points"
   


   ! Real Test
   call random_number(r_points)
   r_points = (r_points(:,:) * (urange + 1 - lrange)) + lrange
   rpsr = (r_points(:,:) * (urange + 1 - lrange)) + lrange
   rps = r_points
   write(0,*) ""
   write(0,*) "=================="
   write(0,*) "Real Test:"
   write(0,*) " Real Shape: ", shape(r_points), " nElms: ", size(r_points)
   write(0,*) ""

   rTree => mpas_kd_construct(r_points)

end subroutine test1

subroutine test2(ndims, lrange, urange)

   implicit none
   integer, intent(in), value :: ndims
   integer, intent(in), value :: lrange, urange

   real, dimension(ndims) :: arry1, arry2

   real :: distance

   call random_number(arry1)
   call random_number(arry2)

   arry1 = (arry1(:) * (urange + 1 - lrange)) + lrange
   arry2 = (arry2(:) * (urange + 1 - lrange)) + lrange
   
   write(0,*) ""
   write(0,*) "arry 1: ", arry1
   write(0,*) ""
   write(0,*) "arry 2: ", arry2
   write(0,*) ""

   distance = sum(arry1(:) - arry2(:))**2

   write(0,*) "Distance: ", distance


end subroutine test2


end program mpas_kd_tester
