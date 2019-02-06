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

      class(*), dimension(:), pointer :: data
   end type kdnode

! Subroutine and function interfaces
   interface mpas_kd_insert
      module procedure  mpas_kd_insert_integer
      module procedure  mpas_kd_insert_real
   end interface mpas_kd_insert

   interface mpas_kd_construct
      module procedure mpas_kd_construct_tree_integer
      module procedure mpas_kd_construct_tree_real
   end interface mpas_kd_construct
   
   interface mpas_kd_search
      module procedure mpas_kd_search_integer
      module procedure mpas_kd_search_real
   end interface mpas_kd_search

   interface mpas_kd_remove
      module procedure mpas_kd_remove_integer
      module procedure mpas_kd_remove_real 
   end interface mpas_kd_remove


! Operator Interfaces
   interface assignment(=)
      module procedure assign_int_to_node
      module procedure assign_real_to_node
      module procedure assign_node_to_node
   end interface assignment(=)

   interface operator (==)
      module procedure mpas_kd_real_equality
      module procedure mpas_kd_integer_equality
   end interface operator (==)

   interface operator (>)
      module procedure mpas_kd_real_gt
      module procedure mpas_kd_integer_gt
   end interface operator (>)

   interface operator (<)
      module procedure mpas_kd_real_lt
      module procedure mpas_kd_integer_lt
   end interface operator (<)

   interface operator (>=)
      module procedure mpas_kd_real_gte
      module procedure mpas_kd_integer_gte
   end interface operator (>=)

   interface operator (<=)
      module procedure mpas_kd_real_lte
      module procedure mpas_kd_integer_lte
   end interface operator (<=)

   interface bubblesort
      module procedure bubbleSort_integer
      module procedure bubbleSort_real
   end interface bubblesort

   contains


! mpas_kd_insert()
! 
   subroutine mpas_kd_insert_integer(kdtree, val)

      implicit none
      ! Input Variables
      type(kdnode), intent(inout), pointer :: kdtree
      integer, dimension(:), intent(in) :: val

   end subroutine mpas_kd_insert_integer

   subroutine mpas_kd_insert_real(kdtree, val)

      implicit none
      ! Input Variables
      type(kdnode), intent(inout), pointer :: kdtree
      real, dimension(:), intent(in) :: val

   end subroutine mpas_kd_insert_real


! mpas_kd_construct()
!
   recursive function mpas_kd_construct_tree_integer(points, dim) result(tree)
      
      implicit none
      ! Input Varaibles
      integer, dimension(:,:), intent(inout) :: points
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

     ! write(0,*) "CONSTRUCT INT: Number of points: ", npoints
     ! write(0,*) "CONSTRUCT INT: Number of dims: ", ndims

      if (npoints < 1) then
         !write(0,*) "No points were passed in"
         tree => null()
         return
      endif

      if ( .NOT. present(dim)) then
         d = 0
      else
         d = dim
      endif

      d = mod(dim, ndims) + 1
      write(0,*) d

      allocate(tree)
      median = (1 + npoints) / 2 
      call bubbleSort(points, d)
      
      tree = points(:,median)

      !write(0,*) " Sorted Points: ", points(d,:)
      !write(0,*) " Median: ", points(d, median)
   
      if (median >= 1) then ! Go left
         !write(0,*) "Went left"
         tree % left = mpas_kd_construct(points(:,1:median-1), d)
      else
         tree % left => null() 
      endif
      
      if (median+1 <= npoints) then ! Go right
         !write(0,*) "Went right"
         tree % right = mpas_kd_construct(points(:,median+1:npoints), d)
      else
         tree % right => null()
      endif

   end function mpas_kd_construct_tree_integer

   recursive function mpas_kd_construct_tree_real(points, dim) result(tree)

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

     ! write(0,*) "CONSTRUCT INT: Number of points: ", npoints
     ! write(0,*) "CONSTRUCT INT: Number of dims: ", ndims

      if(npoints < 1) then
         tree => null()
         return
      endif

      if ( .NOT. present(dim)) then
         d = 0
      else
         d = dim
      endif

      d = mod(dim, ndims) + 1

      allocate(tree)
      median = (1 + npoints) / 2 
      call bubbleSort(points, d)
      
      tree = points(:,median)

      !write(0,*) " Sorted Points: ", points(d,:)
      !write(0,*) " Median: ", points(d, median)
   
      if (median >= 1) then ! Go left
         !write(0,*) "Went left"
         tree % left = mpas_kd_construct(points(:,1:median-1), d)
      else
         tree % left => null() 
      endif
      
      if (median+1 <= npoints) then ! Go right
         !write(0,*) "Went right"
         tree % right = mpas_kd_construct(points(:,median+1:npoints), d)
      else
         tree % right => null()
      endif

   end function mpas_kd_construct_tree_real


! mpas_kd_find()
! Search returns the closet value 
   recursive function mpas_kd_search_integer(kdtree, point, dim) result(result)
      
      implicit none

      ! Input Variables
      type(kdnode), intent(in) :: kdtree
      integer, dimension(:), intent(in) :: point
      integer, value :: dim
      
      ! Return Value
      integer :: result
      integer :: ndims
      integer :: d

      ndims = size(kdtree % data, dim=1)

      d = mod(dim, ndims) + 1


      if(associated(kdtree % right)) then
      endif
      if(associated(kdtree % left)) then

      endif


   end function mpas_kd_search_integer

   function mpas_kd_search_real(kdtree, point, dim) result(result)
      
      implicit none

      ! Input Variables
      type(kdnode), intent(in) :: kdtree
      real, dimension(:), intent(in) :: point
      integer, value :: dim

      ! Return Value
      real :: result
      integer :: ndims
      integer :: d

      ndims = size(kdtree%data, dim=1)

      d = mod(dim, ndims) + 1

   end function mpas_kd_search_real


! mpas_kd_del()
!
   function mpas_kd_remove_integer(kdtree, val) result(ierr)

      implicit none
      ! Input Variables
      type (kdnode), intent(inout) :: kdtree
      integer, dimension(:), intent(in) :: val

      ! Return Value
      integer :: ierr

   end function mpas_kd_remove_integer

   function mpas_kd_remove_real(kdtree, val) result(ierr)

      implicit none
      ! Input Variables
      type (kdnode), intent(inout) :: kdtree
      real, dimension(:), intent(in) :: val

      ! Return Value
      integer :: ierr

   end function mpas_kd_remove_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Operators
!

! = - Assignment
!
   subroutine assign_int_to_node(LHS, RHS)
      implicit none
      ! Input Variables
      type(kdnode), intent(inout) :: LHS
      integer, dimension(:), intent(in) :: RHS

   end subroutine assign_int_to_node
   
   subroutine assign_real_to_node(LHS, RHS)
      implicit none
      type(kdnode), intent(inout) :: LHS
      real, dimension(:), intent(in) :: RHS


   end subroutine assign_real_to_node

   subroutine assign_node_to_node(LHS, RHS)
      implicit none
      type(kdnode), intent(inout) :: LHS
      type(kdnode), intent(in) :: RHS

   end subroutine assign_node_to_node

! == - Equivalence
!
   function mpas_kd_integer_equality(ival, tval) result(result)

      implicit none
      ! Input Variable
      integer, dimension(:), intent(in) :: ival
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_integer_equality

   function mpas_kd_real_equality(rval, tval) result(result)

      implicit none
      ! Input Variable
      real, dimension(:), intent(in) :: rval
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_real_equality


   ! > - Greater Than
   !
   function mpas_kd_integer_gt(ival, tval) result(result)

      implicit none
      ! Input Variable
      integer, dimension(:), intent(in) :: ival
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_integer_gt

   function mpas_kd_real_gt(rval, tval) result(result)

      implicit none
      ! Input Variable
      real, dimension(:), intent(in) :: rval
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_real_gt

! < - Less Than
!
   function mpas_kd_integer_lt(ival, tval) result(result)

      implicit none
      ! Input Variable
      integer, dimension(:), intent(in) :: ival
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_integer_lt

   function mpas_kd_real_lt(rval, tval) result(result)

      implicit none
      ! Input Variable
      real, dimension(:), intent(in) :: rval
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_real_lt

! >= - Greater Than or Equal
!
   function mpas_kd_integer_gte(ival, tval) result(result)

      implicit none
      ! Input Variable
      integer, dimension(:), intent(in) :: ival
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_integer_gte

   function mpas_kd_real_gte(rval, tval) result(result)

      implicit none
      ! Input Variable
      real, dimension(:), intent(in) :: rval
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_real_gte

! <= - Less Than or Equal
!
   function mpas_kd_integer_lte(ival, tval) result(result)

      implicit none
      ! Input Variable
      integer, dimension(:), intent(in) :: ival
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_integer_lte

   function mpas_kd_real_lte(rval, tval) result(result)

      implicit none
      ! Input Variable
      real, dimension(:), intent(in) :: rval
      type(kdnode), intent(in) :: tval

      ! Return Value
      logical :: result

   end function mpas_kd_real_lte

!!!! Possible routines:

   ! mpas_kd_balance()
   ! mpas_kd_nn(num) - Find the num nearest neighbors to a point
   !  num - The number of nearest neighbors desired
   ! mpas_kd_range_search() - Find a range across axis 
   ! Ex: Return the nodes that are between x1, x2 and y1, y2

!!!! Sorts

   subroutine bubbleSort_integer(array, dim)
      implicit none
      ! Input Variable
      integer, dimension(:,:), intent(inout) :: array
      integer, value :: dim

      ! Local Variable
      integer, dimension(size(array, dim=1)) :: swap
      integer :: npoints
      integer :: i, j
      logical :: sorted

      npoints = size(array, dim=2)

      do i = 1, npoints, 1
         do j = 1, npoints - i, 1
            if(array(dim, j + 1) < array(dim, j)) then
               swap = array(:, j)
               array(:, j) = array(:, j + 1)
               array(:, j + 1) = swap
            endif
         enddo
      enddo

   end subroutine bubbleSort_integer

   subroutine bubbleSort_real(array, dim)

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

   end subroutine bubbleSort_real


end module mpas_kd_tree
