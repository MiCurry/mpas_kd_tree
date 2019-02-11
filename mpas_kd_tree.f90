module mpas_kd_tree
    
   implicit none

   private

   public :: kdnode

   ! Public Subroutines
   public :: mpas_kd_insert
   public :: mpas_kd_construct
   public :: mpas_kd_search
   public :: mpas_kd_remove

   ! Public Operators
   public :: operator(==)
   public :: operator(>)
   public :: operator(<)
   public :: operator(>=)
   public :: operator(<=)

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
   subroutine mpas_kd_insert(kdtree, val)

      implicit none
      ! Input Variables
      type(kdnode), intent(inout), pointer :: kdtree
      real, dimension(:), intent(in) :: val

   end subroutine mpas_kd_insert


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpas_kd_construct()
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive function mpas_kd_construct(points, dim) result(tree)

      implicit none
      ! Input Varaibles
      real, dimension(:,:), intent(inout) :: points
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

      d = mod(dim, ndims) + 1

      median = (1 + npoints) / 2 
      call bubbleSort(points, d) ! Sort the points

     ! write(0,*) " Sorted Points: ", points(:,:)
     ! write(0,*) " Median point: ", points(:, median)
     ! write(0,*) " "

      allocate(tree) ! Allocate the node
      allocate(tree % data(ndims)) ! Allocate the data for that node
      tree % data = points(:,median)

      !write(0,*) "We allocated the point: ", points(:, median), tree % data


      !write(0,*) " Median: ", points(d, median)
   
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
   !! mpas_kd_find()
   !! Search returns the closet value 
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

      ! Return Value
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
         min_d = current_distance
         result = kdtree % data(:) 
      endif

      if ( point(d) > kdtree % data(d) ) then

         if ( associated( kdtree % right ) ) then ! Search Right
            call mpas_kd_search(kdtree % right, point, result, min_d, d)
         endif
         if ( sum(point(d) - kdtree % data(:))**2 <= min_d .AND. associated(kdtree % left)) then 
            call mpas_kd_search(kdtree % left, point, result, min_d, d)
         endif

      else if ( point(d) < kdtree % data(d) ) then ! Search Left
         if ( associated( kdtree % left) ) then
            call mpas_kd_search(kdtree % left, point, result, min_d, d)
         endif
         if ( sum(point(:) - kdtree % data(:))**2 <= min_d .AND. associated(kdtree % right )) then
            call mpas_kd_search(kdtree % right, point, result, min_d, d)
         endif

      else ! Go both
         if(associated(kdtree % right)) call mpas_kd_search(kdtree % right, point, result, min_d, d+1)
         if(associated(kdtree % left)) call mpas_kd_search(kdtree % left, point, result, min_d, d+1)
      endif

   end subroutine mpas_kd_search


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! mpas_kd_remove
   !!
   !!
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mpas_kd_remove(kdtree, val) result(ierr)

      implicit none
      ! Input Variables
      type (kdnode), intent(inout) :: kdtree
      real, dimension(:), intent(in) :: val

      ! Return Value
      integer :: ierr

   end function mpas_kd_remove


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


end module mpas_kd_tree
