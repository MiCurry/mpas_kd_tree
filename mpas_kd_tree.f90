module mpas_kd_tree
    
   implicit none

   private

   public :: kdnode

   ! Public Subroutines
   public :: mpas_kd_insert
   public :: mpas_kd_construct
   public :: mpas_kd_search
   public :: mpas_kd_remove
   public :: mpas_kd_free

   type kdnode
      type (kdnode), pointer :: left => null()
      type (kdnode), pointer :: right => null()

      real, dimension(:), pointer :: data
   end type kdnode

   contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpas_kd_insert()
   !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive subroutine mpas_kd_insert(kdtree, val, dim)

      implicit none
      ! Input Variables
      type(kdnode), intent(inout), pointer :: kdtree
      real, dimension(:), intent(in) :: val
      integer, optional, value :: dim

      integer :: ndims
      integer :: d

      if ( .NOT. present(dim)) then
         d = 0
      else
         d = dim
      endif

      ndims = size(val)
      d = mod(d, ndims) + 1

      if ( .NOT. associated(kdtree)) then
         allocate(kdtree)
         allocate(kdtree % data(ndims))
         kdtree % left => null()
         kdtree % right => null()
         kdtree % data = val
         return
      endif

      if ( val(d) > kdtree % data(d)) then
         call mpas_kd_insert(kdtree % right , val, d)
      else
         call mpas_kd_insert(kdtree % left, val, d)
      endif

   end subroutine mpas_kd_insert


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpas_kd_construct()
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive function mpas_kd_construct(points, dim) result(tree)

      implicit none
      ! Input Varaibles
      real, dimension(:,:) :: points
      integer, optional, value :: dim

      ! Return Value
      type (kdnode), pointer :: tree

      ! Local Variables
      integer :: npoints
      integer :: ndims
      integer :: d
      integer :: median

      ndims = size(points, dim=1)
      npoints = size(points, dim=2)

     ! write(0,*)
     ! write(0,*) "CONSTRUCT INT: Number of points: ", npoints
     ! write(0,*) "CONSTRUCT INT: Number of dims: ", ndims

      if(npoints < 1) then
         tree => null()
         return
      endif

      ! If we have no dimensions
      if ( .NOT. present(dim)) then
         d = 0
      else
         d = dim
      endif

      d = mod(d, ndims) + 1

      median = (1 + npoints) / 2 
      
      ! Sort the points
      !call bubbleSort(points, d) ! Sort the points
      call quickSort(points, d, 1, npoints) ! Sort the points

      allocate(tree) ! Allocate the node
      allocate(tree % data(ndims)) ! Allocate the data for that node
      tree % data = points(:,median)

      if (median >= 1) then ! Go left
         tree % left => mpas_kd_construct(points(:,1:median-1), d)
      else
         tree % left => null() 
      endif
      
      if (median+1 <= npoints) then ! Go right
         tree % right => mpas_kd_construct(points(:,median+1:npoints), d)
      else
         tree % right => null()
      endif

   end function mpas_kd_construct


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpas_kd_find(kdtree, point, result, min_d, dim)
   !! - Find the closet value to point(:)
   !!
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive subroutine mpas_kd_search(kdtree, point, result, min_d, dim)
      
      implicit none

      ! Input Variables
      type(kdnode), pointer, intent(in) :: kdtree
      real, dimension(:), intent(in) :: point
      real, dimension(:), intent(out) :: result
      real, intent(inout) :: min_d
      integer, optional, value :: dim 

      ! Local Values
      integer :: ndims
      integer :: d
      real :: best_distance, current_distance

      ! Check to see if result and point have the at least the same dimensions?
      ndims = size(kdtree % data, dim=1)

      if ( .not. present(dim)) then
         d = 0
      else
         d = dim
      endif

      d = mod(d, ndims) + 1

      current_distance = sum( kdtree % data(:) - point(:) )**2
      if ( current_distance < min_d ) then
         !write(0,*) "min: ", min_d, "current_dist: ", current_distance
         min_d = current_distance
         result = kdtree % data(:) 
         !write(0,*) "Result: ", result
      endif

      if ( point(d) > kdtree % data(d) ) then
         if ( associated( kdtree % right ) ) then ! Search Right
            call mpas_kd_search(kdtree % right, point, result, min_d, d)
         endif
         if ( sum(point(:) - kdtree % data(:))**2 <= min_d .AND. associated(kdtree % left)) then 
            call mpas_kd_search(kdtree % left, point, result, min_d, d)
         endif

      else if ( point(d) < kdtree % data(d) ) then ! Search Left
         if ( associated( kdtree % left) ) then
            call mpas_kd_search(kdtree % left, point, result, min_d, d)
         endif
         if ( sum(point(:) - kdtree % data(:))**2 <= min_d .AND. associated(kdtree % right)) then
            call mpas_kd_search(kdtree % right, point, result, min_d, d)
         endif

      else ! Go both
         if(associated(kdtree % right)) call mpas_kd_search(kdtree % right, point, result, min_d, d)
         if(associated(kdtree % left)) call mpas_kd_search(kdtree % left, point, result, min_d, d)
      endif

   end subroutine mpas_kd_search


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpas_kd_remove
   !!
   !!
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpas_kd_remove(kdtree, point) result(ierr)

      implicit none
      ! Input Variables
      type (kdnode), intent(inout) :: kdtree
      real, dimension(:), intent(in) :: point

      ! Return Value
      integer :: ierr

      

   end function mpas_kd_remove

   recursive subroutine mpas_kd_free(kdtree)

      implicit none
      type(kdnode), pointer :: kdtree


      if (associated(kdtree % left)) then
         call mpas_kd_free(kdtree % left)
      endif

      if (associated(kdtree % right)) then
         call mpas_kd_free(kdtree % right)
      endif

      deallocate(kdtree % data)
      deallocate(kdtree)

   end subroutine mpas_kd_free


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sorts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine bubbleSort(array, dim)

      implicit none
      ! Input Variable
      real, dimension(:,:), intent(inout) :: array
      integer, value :: dim

      ! Local Variable
      real, dimension(size(array, dim=1)) :: swap
      integer :: npoints
      integer :: i, j, k

      npoints = size(array, dim=2)

      do i = 1, npoints, 1
         k = i
         do j = 1, npoints - k, 1
            if(array(dim, j + 1) < array(dim, j)) then
               swap = array(:, j)
               array(:, j) = array(:, j + 1)
               array(:, j + 1) = swap
            endif
         enddo
      enddo

   end subroutine bubbleSort

   recursive subroutine quickSort(array, dim, arrayStart, arrayEnd)

      implicit none
      ! Input Variable
      real, dimension(:,:) :: array
      integer, intent(in), value :: dim
      integer, intent(in), value :: arrayStart, arrayEnd

      ! Local Variables 
      integer :: ndims, npoints
      real, dimension(size(array, dim=1)) :: temp
      real, dimension(size(array, dim=1)) :: pivot_value

      integer :: l, r, pivot, s


      ndims = size(array, dim=1)
      npoints = arrayEnd

      if ( (arrayEnd - arrayStart) < 1 ) then
         return
      endif

      ! Create the left, right, and start pointers
      l = arrayStart
      r = arrayEnd - 1
      s = l

      pivot = (l+r)/2
      pivot_value = array(:, pivot)

      ! Move the pivot to the far right
      temp(:) = array(:,pivot)
      array(:,pivot) = array(:,arrayEnd)
      array(:,arrayEnd) = temp(:)

      do while ( .TRUE. )
         ! Advance the left pointer until it is a value less then our pivot_value(dim)
         do while ( array(dim, l) < pivot_value(dim) )
            l = l + 1
         enddo

         ! Advance the right pointer until it is a value more then our pivot_value(dim)
         do while ( r > 0 .AND. array(dim, r) > pivot_value(dim) )
            r = r - 1
         enddo

         ! l >= r so return and quicksort the left and right subtrees
         if ( l >= r ) then 
            exit
         else ! Swap elements about the pivot
            temp = array(:,l)
            array(:,l) = array(:,r)
            array(:,r) = temp
         endif
      enddo

      ! Move the pivot to where the l ended up
      temp(:) = array(:,l)
      array(:,l) = array(:,arrayEnd)
      array(:,arrayEnd) = temp(:)

      !Quick Sort on the lower partition
      call quickSort(array(:,:), dim, s, l-1)
      
      !Quick sort on the upper partition
      call quickSort(array(:,:), dim, l+1, arrayEnd)

   end subroutine quicksort

   subroutine insertionSort(array, dim, startArray, endArray)

      implicit none

      real, dimension(:,:), intent(inout) :: array
      integer, intent(in), value :: dim
      integer, intent(in), value :: startArray
      integer, intent(in), value :: endArray

      real, dimension(dim) :: temp
      integer :: i, j

      do i = startArray, endArray, 1
         do j = startArray, endArray, 1
            if ( array(dim, i) > array(dim, j) ) then
               temp(:) = array(:,i)
               array(:,i) = array(:,j)
               array(:,j) = temp(:)
            endif
         enddo
      enddo

   end subroutine insertionSort

end module mpas_kd_tree
