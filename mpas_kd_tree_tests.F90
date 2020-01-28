program mpas_kd_tester

   use mpas_kd_tree, only : mpas_kd_type
   use mpas_kd_tree, only : mpas_kd_construct, mpas_kd_search
   use mpas_kd_tree, only : quicksort
   use mpas_kd_tree, only : mpas_kd_free

   implicit none

#ifdef SINGLE_PRECISION
   integer, parameter :: RKIND = selected_real_kind(6)
#else
   integer, parameter :: RKIND = selected_real_kind(12)
#endif

   integer :: lrange
   integer :: urange

   lrange = -100000
   urange = 100000

   call test1()
   call test2(5000, 1)
   call test2(5000, 3)
   call test2(5000, 3)
   call test2(5000, 7)
   call quicksort_test(5000, 1, 1)
   call quicksort_test(5000, 3, 1)
   call quicksort_test(5000, 3, 2)
   call quicksort_test(5000, 3, 3)

contains

! Test create, search and free 1, 2, 3 and 7 elements
subroutine test1()

   implicit none

   type(mpas_kd_type), dimension(:), pointer :: nodes
   type(mpas_kd_type), pointer :: tree => null()
   type(mpas_kd_type), pointer :: res => null()

   write(0,*) "MPAS_KD_TREE Tests Started"

   !
   ! Trying to free an unassociated tree
   !
   call mpas_kd_free(tree)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
   endif

   !
   ! Testing mpas_kd_construct with 1 node
   !
   write(0,'(a)',advance='no') "Testing mpas_kd_construct with 1 node - "

   allocate(nodes(1))
   allocate(nodes(1) % point(3))
   nodes(1) % point(:) = (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/)
   nodes(1) % id = 1

   tree => mpas_kd_construct(nodes, 3)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
   endif

   ! Root node should contain our point
   if (all(tree % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Root node did not contain the correct point"
      stop
   endif

   ! Check that the leaf nodes are null()
   if (associated(tree % left) .or. associated(tree % right)) then
       write(0,*) "ERROR: Left and right nodes were associated when they should not have been"
       stop
   endif


   call mpas_kd_search(tree, (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/), res)
   if (all(res % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % id /= 1) then
      write(0,*) "ERROR: Cell did not match the expected id"
   endif

   call mpas_kd_search(tree, (/1000.0_RKIND, 1000.0_RKIND, 1000.0_RKIND/), res)
   if (all(res % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % id /= 1) then
      write(0,*) "ERROR: Cell did not match the expected id"
   endif

   write(0,*) ""
   write(0,*) "Testing Free: "
   call mpas_kd_free(tree)
   if (associated(tree)) then
      write(0,*) "Tree was associated"
      stop 1
   endif
   write(0,*) "Tree was not associated"
   deallocate(nodes)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
   endif

   write(0,*) "PASS"

   ! 
   ! Testing mpas_kd_construct with 2 nodes
   !
   write(0,'(a)',advance='no') "Testing mpas_kd_construct with 2 nodes - "
   allocate(nodes(2))

   allocate(nodes(1) % point(3))
   nodes(1) % point(:) = (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/)
   nodes(1) % id = 1

   allocate(nodes(2) % point(3))
   nodes(2) % point(:) = (/-1.0_RKIND, -1.0_RKIND, -1.0_RKIND/)
   nodes(2) % id = 2

   tree => mpas_kd_construct(nodes, 3)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
   endif

   call mpas_kd_search(tree, (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/), res)
   if (all(res % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % id /= 1) then
      write(0,*) "ERROR: Cell did not match the expected id"
   endif

   call mpas_kd_search(tree, (/-1.0_RKIND, -1.0_RKIND, -1.0_RKIND/), res)
   if (all(res % point(:) /= (/-1.0_RKIND, -1.0_RKIND, -1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct"
   endif
   if (res % id /= 2) then
      write(0,*) "ERROR: Cell did not match the expected id"
   endif

   call mpas_kd_free(tree)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
   endif

   deallocate(nodes)

   write(0,*) "PASS"

   !
   ! Testing mpas_kd_construct with 3 nodes
   !
   write(0,'(a)',advance='no') "Testing mpas_kd_construct with 3 nodes - "
   allocate(nodes(3))

   allocate(nodes(1) % point(3))
   nodes(1) % point(:) = (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/)
   nodes(1) % id = 1

   allocate(nodes(2) % point(3))
   nodes(2) % point(:) = (/2.0_RKIND, 2.0_RKIND, 2.0_RKIND/)
   nodes(2) % id = 2

   allocate(nodes(3) % point(3))
   nodes(3) % point(:) = (/3.0_RKIND, 3.0_RKIND, 3.0_RKIND/)
   nodes(3) % id = 3

   tree => mpas_kd_construct(nodes, 3)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
   endif

   ! Testing if the root nodes 2nd row's children are null()
   if (associated(tree % left % left) .or. associated(tree % left % right)) then
      write(0,*) "ERROR: One of the left child's leaf nodes was associated when it should have not been"
      stop 1
   endif

   if (associated(tree % right % left) .or. associated(tree % right % right)) then
      write(0,*) "ERROR: One of the right child's leaf nodes was associated when it should have not been"
      stop 1
   endif

   call mpas_kd_search(tree, (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/), res)
   if (all(res % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct (1.0, 1.0, 1.0)"
      stop 1
   endif
   if (res % id /= 1) then
      write(0,*) "ERROR: Cell did not match the expected id 1"
      stop 1
   endif

   call mpas_kd_search(tree, (/2.0_RKIND, 2.0_RKIND, 2.0_RKIND/), res)
   if (all(res % point(:) /= (/2.0_RKIND, 2.0_RKIND, 2.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct (2.0, 2.0, 2.0)", res % point(:)
      stop 1
   endif
   if (res % id /= 2) then
      write(0,*) "ERROR: Cell did not match the expected id 2"
      stop 1
   endif

   call mpas_kd_search(tree, (/3.0_RKIND, 3.0_RKIND, 3.0_RKIND/), res)
   if (all(res % point(:) /= (/3.0_RKIND, 3.0_RKIND, 3.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct (3.0, 3.0, 3.0)"
      stop 1
   endif
   if (res % id /= 3) then
      write(0,*) "ERROR: Cell did not match the expected id 3"
      stop 1
   endif

   call mpas_kd_free(tree)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
      stop 1
   endif

   deallocate(nodes)

   write(0,*) "PASS"

   !
   ! Testing mpas_kd_construct with 1 node - again
   !
   write(0,'(a)',advance='no') "Testing mpas_kd_construct with 1 node (again) - "

   allocate(nodes(1))
   allocate(nodes(1) % point(3))
   nodes(1) % point(:) = (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/)
   nodes(1) % id = 1

   tree => mpas_kd_construct(nodes, 3)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
      stop 1
   endif

   ! Root node should contain our point
   if (all(tree % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Root node did not contain the correct point"
      stop 1
   endif

   ! Check that the leaf nodes are null()
   if (associated(tree % left) .or. associated(tree % right)) then
       write(0,*) "ERROR: Left and right nodes were associated when they should not have been"
      stop 1
   endif


   call mpas_kd_search(tree, (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/), res)
   if (all(res % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct"
      stop 1
   endif
   if (res % id /= 1) then
      write(0,*) "ERROR: Cell did not match the expected id"
      stop 1
   endif

   call mpas_kd_search(tree, (/1000.0_RKIND, 1000.0_RKIND, 1000.0_RKIND/), res)
   if (all(res % point(:) /= (/1.0_RKIND, 1.0_RKIND, 1.0_RKIND/))) then
      write(0,*) "ERROR: Point from search was not correct"
      stop 1
   endif
   if (res % id /= 1) then
      write(0,*) "ERROR: Cell did not match the expected id"
      stop 1
   endif

   call mpas_kd_free(tree)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
      stop 1
   endif

   deallocate(nodes)

   write(0,*) "PASS"

end subroutine test1

subroutine test2(n, ndims)
   
   implicit none

   integer, intent(in), value :: n
   integer, intent(in), value :: ndims

   real(kind=RKIND), dimension(:,:), pointer :: array1 => null()
   real(kind=RKIND), dimension(:,:), pointer :: array2 => null()
   real(kind=RKIND), dimension(ndims) :: min_point
   real(kind=RKIND) :: min_distance, distance, kd_distance

   type(mpas_kd_type), pointer :: tree => null(), res => null()
   type(mpas_kd_type), dimension(:), pointer :: nodes => null()
   integer :: i, j, num_missed = 0

   write(0,*) "Testing with large amounts of points ", n, ndims
   
   allocate(array1(ndims, n))
   allocate(array2(ndims, n))

   ! Create two arrays, one for our tree and one to check
   call random_number(array1(:,:))
   call random_number(array2(:,:))
   array1 =  (array1(:,:) * (urange + 1 - lrange)) + lrange
   array2 =  (array2(:,:) * (urange + 1 - lrange)) + lrange

   allocate(nodes(n))
   do i = 1, n
      allocate(nodes(i) % point(ndims))
      nodes(i) % point(:) = array1(:,i)
      nodes(i) % id = i
   enddo

   tree => mpas_kd_construct(nodes, ndims)
   if (.not. associated(tree)) then
      write(0,*) "ERROR: Tree was unassociated when it should have been associated"
      stop 1
   endif

   write(0,*) "Searching for points that we used to create the tree"
   write(0,*) "Brute forcing answer to confirm correct kd search"
   write(0,*) "This will take a while with large n and higher dimensions..."
   do i = 1, n
      call mpas_kd_search(tree, array1(:,i), res)
      if (all(res % point(:) /= array1(:, i))) then
         write(0,*) "ERROR: Could not find the correct point!"
         stop 1
      endif
      if (res % id /= i) then
         write(0,*) "ERROR: Cell did not equal what it should have!"
         write(0,*) "  Returned Node: ", res % point(:), res % id
         write(0,*) "  Desired id: ", array1(:,i), i
      endif
     
      ! Brute force to see if we get the right answer
      min_distance = huge(min_distance)
      do j = 1, n
         distance = sum((array1(:,j) - array1(:,i))**2)
         if (distance < min_distance) then
            min_point = array1(:,j)
            min_distance = distance
         endif
      enddo

      if (all(min_point(:) /= res % point(:))) then
         write(0,*) "ERROR: The point found via brute force was not the same as the one from the kd search"
         write(0,*) "ERROR: Brute Force Point: ", min_point(:), " KD search point: ", res % point(:)
         stop 1
      endif
   enddo

   write(0,*) "Searching for points that are not in the tree"
   write(0,*) "Brute forcing answer to confirm correct kd search"
   write(0,*) "This will take a while with large n and higher dimensions..."
   do i = 1, n
      call mpas_kd_search(tree, array2(:,i), res, kd_distance)
     
      ! Brute force to see if we get the right answer
      min_distance = huge(min_distance)
      do j = 1, n
         distance = sum((array1(:,j) - array2(:,i))**2)
         if (distance < min_distance) then
            min_point = array1(:,j)
            min_distance = distance
         endif
      enddo

      if (all(min_point(:) /= res % point(:))) then
         write(0,*) "ERROR: The point found via brute force was not the same as the one from the kd search"
         write(0,*) "ERROR: Searched for: ", array2(:,i)
         write(0,*) "ERROR: Brute Force Point: ", min_point(:), "Distance:", min_distance
         write(0,*) "ERROR: KD search point: ", res % point(:), "Distance:", kd_distance
         num_missed = num_missed + 1
      endif
   enddo

   deallocate(array1)
   deallocate(array2)
   
   call mpas_kd_free(tree)
   if (associated(tree)) then
      write(0,*) "ERROR: Tree was associated after calling mpas_kd_free"
      stop 1
   endif

   deallocate(nodes)

   write(0,*) "" 
   write(0,*) "Number of Points Missed: ", num_missed

   write(0,*) "PASS"

end subroutine test2

subroutine quicksort_test(n, ndims, dim)

   implicit none

   integer, intent(in), value :: n
   integer, intent(in), value :: ndims
   integer, intent(in), value :: dim
   type(mpas_kd_type), dimension(:), pointer :: nodes => null()
   real(kind=RKIND), dimension(:,:), pointer :: array => null()
   integer :: i

   write(0,*) ""
   write(0,*) "Starting quicksort tests with:"
   write(0,*) "n: ", n 
   write(0,*) "ndims:", ndims
   write(0,*) "Sort on dim: ", dim


   allocate(nodes(n))
   allocate(array(n,ndims))

   call random_number(array(:,:))
   array(:,:) =  (array(:,:) * (urange + 1 - lrange)) + lrange

   do i = 1, n
      allocate(nodes(i) % point(ndims))
      nodes(i) % point(:) = array(i,:)
   enddo

   call quicksort(nodes, dim, 1, size(nodes), ndims)

   do i = 2, n
      if (.not. (nodes(i - 1) % point(dim) <= nodes(i) % point(dim))) then
         write(0,*) "FAILED! quicksorted returned an unsorted list!"
         stop 1
      endif
   enddo

   do i = 1, n
      deallocate(nodes(i) % point)
   enddo

   deallocate(nodes)
   deallocate(array)

   write(0,*) "Quicksort: PASS"

end subroutine quicksort_test

end program mpas_kd_tester
