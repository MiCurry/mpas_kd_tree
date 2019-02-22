program mpas_kd_tester

   use mpas_kd_tree, only : kdnode
   use mpas_kd_tree, only : mpas_kd_insert, mpas_kd_construct, mpas_kd_search, mpas_kd_remove
   use mpas_kd_tree, only : mpas_kd_free
   use mpas_kd_tree, only : mpas_kd_find_min
   use getoptf

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
            call test3(npoints, ndims, lrange, urange)
         case ('4')
            write(0,*) "Launching test 4"
            call test4(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('5')
            write(0,*) "Launching test 5"
            call test5(ntests, npoints, ndims, lrange, urange)
         case ('6')
            write(0,*) "Launching test 6"
            call test6(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
         case ('m')
            write(0,*) "Launching Minimum Test"
            call test_min(ntests, npoints, ndims, lrange, urange)
         case ('r')
            write(0,*) "Launching remove Test"
            call test_remove(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)
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
   write(0,*) "Num Points: ", npoints
   write(0,*) "Num Dims: ", ndims

   call random_number(arry1)
   call random_number(arry2)

   arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange
   

   !write(0,*) "arry 1: ", arry1
   !write(0,*) ""
   !write(0,*) "arry 2: ", arry2


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

   call mpas_kd_free(rTree)

end subroutine test2

subroutine test3(npoints, ndims, lrange, urange)

   use iso_c_binding, only : c_int

   implicit none

   interface
       subroutine timer_start(timer_id) bind(C)
          use iso_c_binding, only : c_int
          integer (c_int), intent(in), value :: timer_id
       end subroutine timer_start

       subroutine timer_stop(timer_id, sec, nsec) bind(C)
          use iso_c_binding, only : c_int
          integer (c_int), intent(in), value :: timer_id
          integer (c_int), intent(out) :: sec, nsec
       end subroutine timer_stop
   end interface
 
   integer (c_int) :: timer_id, sec, nsec 
   integer, intent(in), value :: npoints, ndims
   integer, intent(in), value :: lrange, urange

   type(kdnode), pointer :: rTree => null()
   real, dimension(ndims, npoints) :: arry1, arry2
   real, dimension(ndims) :: result
   real :: min_d

   integer :: i, r, k, n

   min_d = huge(min_d)

   write(0,*) "Timing test. Size: ", npoints, " dims: ", ndims
   call random_number(arry1)
   call random_number(arry2)

   arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

   write(0,*) "Constructing tree..." 
   call timer_start(0)
   rTree => mpas_kd_construct(arry1) 
   call timer_stop(0, sec, nsec)
   write(0,*) "Tree constructed... time: ", sec, nsec

   
   write(0,*) "Searching for points that are within the tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call mpas_kd_search(rTree, arry1(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Searching for points that are NOT within the tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call mpas_kd_search(rTree, arry2(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Freeing rTree"
   call timer_start(0)
   call mpas_kd_free(rTree)
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "=============================================="
   write(0,*) "Testing timing for inserting nodes"
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call mpas_kd_insert(rTree, arry1(:, i))
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "Tree insertion creation time: ", sec, nsec

   write(0,*) "Searching for points that are within the INSERTED tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call mpas_kd_search(rTree, arry1(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Searching for points that are NOT within the INSERTED tree..." 
   call timer_start(0)
   do i = 1, size(arry1, dim=2)
      call mpas_kd_search(rTree, arry2(:, i), result, min_d)
   enddo
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   write(0,*) "Freeing rTree"
   call timer_start(0)
   call mpas_kd_free(rTree)
   call timer_stop(0, sec, nsec)
   write(0,*) "time: ", sec, nsec

   
end subroutine test3

function search_tree(tree, point)

   implicit none

   type(kdnode), pointer, intent(in) :: tree
   real, dimension(:), intent(in) :: point
   real, dimension(size(point, dim=1)) :: result

   logical search_tree
   real :: min_d

   search_tree = .TRUE.

   min_d = huge(min_d)
   !write(0,*) "New Search", point
   call mpas_kd_search(tree, point, result, min_d)
   if (.NOT. (all(result(:) == point(:)) .AND. min_d == 0.0)) then
      write(0,*) point(:), " was not in the tree, but its cloest point was: ", result, "and the dist: ", sqrt(min_d)
      search_tree = .FALSE.
      return
   endif

end function search_tree

subroutine test4(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   real, dimension(:,:), pointer :: arry1, arry2
   type(kdnode), pointer :: tree => null()

   real :: r_ndims, r_npoints
   integer :: ndims, npoints

   integer :: i, j

   write(0,*) "Running Random tests"
   write(0,*) " running ", ntests, " tests"
   write(0,*) ""

   do i = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) "Test number: ", i
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      tree => mpas_kd_construct(arry1)
      
      !write(0,*) "Arry1: ", arry1

      do j = 1, size(arry1, dim=2), 1
         if( .NOT. search_tree(tree, arry1(:,j))) then
            write(0,*) "This point was not found"
            stop
         endif
      enddo

      deallocate(arry1)
      deallocate(arry2)
      call mpas_kd_free(tree)

   enddo

end subroutine test4

subroutine test5(ntests, npoints, ndims, lrange, urange)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: npoints
   integer, intent(in), value :: ndims
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   integer :: i, j
   real :: r_ndims, r_npoints
   integer :: ndims2, npoints2

   allocate(arry1(ndims, npoints))
   allocate(arry2(ndims, npoints))

   call random_number(arry1(:,:))
   call random_number(arry2(:,:))

   arry1(:,:) = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2(:,:) = (arry2(:,:) * (urange + 1 - lrange)) + lrange

   write(0,*) "Number of points: ", npoints, " Number of dims: ", ndims

   do i = 1, size(arry1, dim=2), 1
      call mpas_kd_insert(tree, arry1(:,i))
   enddo

   write(0,*) "All ", npoints, "were added into the tree!"

   do i = 1, size(arry1, dim=2), 1
      if( .NOT. search_tree(tree, arry1(:,i))) then
         write(0,*) "This point was not found! But it should have been!!"
         stop
      endif
   enddo

   write(0,*) "All points that were inserted to make the tree were succesfully found within the tree!!"
   call mpas_kd_free(tree)

end subroutine test5

subroutine test6(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   real, dimension(:,:), pointer :: arry1, arry2
   type(kdnode), pointer :: tree => null()

   real :: r_ndims, r_npoints
   integer :: ndims, npoints

   integer :: i, j

   write(0,*) "Running Random tests with insert"
   write(0,*) " running ", ntests, " tests"
   write(0,*) ""

   do i = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) "Test number: ", i
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      do j = 1, npoints, 1
         call mpas_kd_insert(tree, arry1(:, j))
      enddo
      
      !write(0,*) "Arry1: ", arry1

      do j = 1, size(arry1, dim=2), 1
         if( .NOT. search_tree(tree, arry1(:,j))) then
            write(0,*) "This point was not found"
            stop
         endif
      enddo

      deallocate(arry1)
      deallocate(arry2)
      call mpas_kd_free(tree)

   enddo

end subroutine test6

subroutine test_min(ntests, npoints, ndims, lrange, urange)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: npoints
   integer, intent(in), value :: ndims
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   real, dimension(2,7) :: arry3
   real, dimension(ndims) :: min_point 
   integer :: i, j
   real :: r_ndims, r_npoints
   real, dimension(ndims) :: minimum
   integer :: ndims2, npoints2

   minimum = huge(minimum)

   arry3 = reshape((/15, 33, 41, 66, 35, 92, 46, 53, 80, 65, 50, 76, 75, 86, 81, 98/), shape(arry3))

   allocate(arry1(ndims, npoints))
   allocate(arry2(ndims, npoints))

   call random_number(arry1(:,:))
   call random_number(arry2(:,:))

   arry1(:,:) = (arry1(:,:) * (urange + 1 - lrange)) + lrange
   arry2(:,:) = (arry2(:,:) * (urange + 1 - lrange)) + lrange

   tree => mpas_kd_construct(arry1(:,:))
   
   call mpas_kd_find_min(tree, min_point, 1)

   write(0,*) "Found min along dimension 1 : ", min_point(1)
   write(0,*) "Minimum along dimension 1: ", minval(arry1(1,:))


end subroutine test_min

recursive subroutine print_tree(tree, left_or_right, level)

   implicit none
   type(kdnode), intent(inout), pointer :: tree
   integer, value :: left_or_right
   integer, value :: level

   if ( .NOT. associated(tree) .OR. .NOT. associated(tree % data)) then
      write(0,*) " I AM AN UNASSOCIATED NODE"
      return
   endif

   if (associated(tree) .AND. associated(tree % data)) then
      if (left_or_right == -1) then
         write(0,*) "Level: ", level, "Left: ", tree % data(:) 
      elseif (left_or_right == 0) then
         write(0,*) "Level: ", level, "Root: ", tree % data(:)
      elseif (left_or_right == 1) then
         write(0,*) "Level: ", level, "Right: ", tree % data(:)
      endif
   endif

   if ( associated(tree % left)) call print_tree(tree % left, -1, level + 1)
   if ( associated(tree % right)) call print_tree(tree % right, 1, level + 1)

end subroutine print_tree

subroutine test_remove(ntests, lrange, urange, npoints_l, npoints_u, ndims_l, ndims_u)

   implicit none

   integer, intent(in), value :: ntests
   integer, intent(in), value :: lrange
   integer, intent(in), value :: urange
   integer, intent(in), value :: npoints_l
   integer, intent(in), value :: npoints_u
   integer, intent(in), value :: ndims_l
   integer, intent(in), value :: ndims_u

   type(kdnode), pointer :: tree => null()
   real, dimension(:,:), pointer :: arry1, arry2
   real, dimension(ndims) :: min_point 
   real, dimension(2,7) :: arry3
   integer :: i, j
   real :: r_ndims, r_npoints
   real, dimension(ndims) :: minimum
   integer :: ndims2, npoints2

   do j = 1, ntests, 1
      call random_number(r_ndims)
      call random_number(r_npoints)
      ndims = int((r_ndims * (ndims_u + 1 - ndims_l)) + ndims_l)
      npoints = int((r_npoints * (npoints_u + 1 - npoints_l)) + npoints_l)
      write(0,*) "Test number: ", j
      write(0,*) "ndims: ", ndims, " npoints: ", npoints

      allocate(arry1(ndims, npoints))
      allocate(arry2(ndims, npoints))

      call random_number(arry1(:,:))
      call random_number(arry2(:,:))
      arry1 = (arry1(:,:) * (urange + 1 - lrange)) + lrange
      arry2 = (arry2(:,:) * (urange + 1 - lrange)) + lrange

      tree => mpas_kd_construct(arry1(:,:))

      do i = 1, size(arry1, dim=2), 1
         if(associated(tree)) then
            if(.NOT. mpas_kd_remove(tree, arry1(:,i))) then
               write(0,*) "We were not able to remove the point: ", arry1(:,i)
               stop
            else
            endif
         else
            write(0,*) "The tree is no long associated"
         endif
      enddo
   enddo

end subroutine

end program mpas_kd_tester
