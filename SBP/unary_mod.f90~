module unary_mod 

  use precision_vars

  private 
  public  :: amubdg , aplbdg, amub, aplb, aplb1, csort, qsorti, qsortd, transp,  &
            & csrcsc, clncsr, csort_block, clean_CSR,                            &
            & a_bl_mu_b_bl, a_sc_mu_b_bl, a_bl_mu_b_sc, a_bl_pl_b_bl1, a_bl_mu_x,&
            & a_bl_pl_b_bl_scal, a_bl_pl_b_bl_diag

contains

  !                          S P A R S K I T                             c
  !----------------------------------------------------------------------c
  !                     UNARY SUBROUTINES MODULE                         c
  !----------------------------------------------------------------------c
  ! contents:                                                            c
  !----------                                                            c
  ! submat : extracts a submatrix from a sparse matrix.                  c
  ! filter : filters elements from a matrix according to their magnitude.c
  ! filterm: same as above, but for the MSR format                       c
  ! csort  : sorts the elements in increasing order of columns           c
  ! clncsr : clean up the CSR format matrix, remove duplicate entry, etc c
  ! transp : in-place transposition routine (see also csrcsc in formats) c
  ! copmat : copy of a matrix into another matrix (both stored csr)      c
  ! msrcop : copies a matrix in MSR format into a matrix in MSR format   c
  ! getelm : returns a(i,j) for any (i,j) from a CSR-stored matrix.      c
  ! getdia : extracts a specified diagonal from a matrix.                c
  ! getl   : extracts lower triangular part                              c
  ! getu   : extracts upper triangular part                              c
  ! levels : gets the level scheduling structure for lower triangular    c
  !          matrices.                                                   c
  ! amask  : extracts     C = A mask M                                   c
  ! rperm  : permutes the rows of a matrix (B = P A)                     c
  ! cperm  : permutes the columns of a matrix (B = A Q)                  c
  ! dperm  : permutes both the rows and columns of a matrix (B = P A Q ) c
  ! dperm1 : general extractiob routine (extracts arbitrary rows)        c
  ! dperm2 : general submatrix permutation/extraction routine            c
  ! dmperm : symmetric permutation of row and column (B=PAP') in MSR fmt c
  ! dvperm : permutes a real vector (in-place)                           c
  ! ivperm : permutes an integer vector (in-place)                       c
  ! retmx  : returns the max absolute value in each row of the matrix    c
  ! diapos : returns the positions of the diagonal elements in A.        c
  ! extbdg : extracts the main diagonal blocks of a matrix.              c
  ! getbwd : returns the bandwidth information on a matrix.              c
  ! blkfnd : finds the block-size of a matrix.                           c
  ! blkchk : checks whether a given integer is the block size of A.      c
  ! infdia : obtains information on the diagonals of A.                  c
  ! amubdg : gets number of nonzeros in each row of A*B (as well as NNZ) c
  ! aplbdg : gets number of nonzeros in each row of A+B (as well as NNZ) c
  ! rnrms  : computes the norms of the rows of A                         c
  ! cnrms  : computes the norms of the columns of A                      c
  ! roscal : scales the rows of a matrix by their norms.                 c
  ! coscal : scales the columns of a matrix by their norms.              c
  ! addblk : Adds a matrix B into a block of A.                          c
  ! get1up : Collects the first elements of each row of the upper        c
  !          triangular portion of the matrix.                           c
  ! xtrows : extracts given rows from a matrix in CSR format.            c
  ! csrkvstr:  Finds block row partitioning of matrix in CSR format      c
  ! csrkvstc:  Finds block column partitioning of matrix in CSR format   c
  ! kvstmerge: Merges block partitionings, for conformal row/col pattern c
  ! qsortd : Quicksort:REAL    algrm. Sorts number array. Accending orderc
  ! qsorti : Quicksort:Integer algrm. Sorts number array. Accending orderc
  !----------------------------------------------------------------------c

  !=============================================================================80

  !----------------------------------------------------------------------c
  !                          S P A R S K I T                             c
  !----------------------------------------------------------------------c
  !        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c
  !----------------------------------------------------------------------c
  ! amub   :   computes     C = A*B                                      c
  ! aplb   :   computes     C = A+B                                      c
  ! aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
  ! aplsb  :   computes     C = A + s B                                  c
  ! aplsb1 :   computes     C = A+sB  [Sorted version: A, B, C sorted]   c
  ! apmbt  :   Computes     C = A +/- transp(B)                          c
  ! aplsbt :   Computes     C = A + s * transp(B)                        c
  ! diamua :   Computes     C = Diag * A                                 c
  ! amudia :   Computes     C = A* Diag                                  c
  ! aplsca :   Computes     A:= A + s I    (s = scalar)                  c
  ! apldia :   Computes     C = A + Diag.                                c
  !----------------------------------------------------------------------c
  !----------------------------------------------------------------------c

  !=============================================================================80
  !                 S P A R S K I T _ B L O C K                          c
  !----------------------------------------------------------------------c
  !        BASIC LINEAR ALGEBRA FOR BLOCK SPARSE MATRICES.               c
  !----------------------------------------------------------------------c
  !                                                                      c
  ! amub   :   computes     C = A*B                                      c
  ! aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
  ! diamua :   Computes     C = Diag * A                                 c
  ! apldia :   Computes     C = A + Diag.                                c
  ! aplsca :   Computes     A:= A + s I    (s = scalar)                  c
  ! amub   :   computes     C = sc(A)*B                                  c


  subroutine submat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao) 
  
    integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*) 
    real(wp) a(*),ao(*) 
    !-----------------------------------------------------------------------
    ! This routine extracts the submatrix A(i1:i2,j1:j2) and puts the result in          
    ! matrix ao,iao,jao
    !
    ! In place: ao,jao,iao may be the same as a,ja,ia.                  
    !
    !--------------                                                         
    ! On input                                                              
    !---------                                                              
    ! n     = row dimension of the matrix                                
    ! i1,i2 = two integers with i2  >=  i1 indicating the range of rows to b
    !          extracted.                                                   
    ! j1,j2 = two integers with j2  >=  j1 indicating the range of columns  
    !         to be extracted.                                              
    !         * There is no checking whether the input values for i1, i2, j1
    !           j2 are between 1 and n.                                     
    ! a,                                                                    
    ! ja,                                                                   
    ! ia    = matrix in compressed sparse row format.                       
    !                                                                       
    ! job   = job indicator: if job  /=  1 then the real values in a ar
    !         extracted, only the column indices (i.e. data structure) are. 
    !         otherwise values as well as column indices are extracted...   
    !                                                                       
    ! On output                                                             
    !--------------                                                         
    ! nr      = number of rows of submatrix                               
    ! nc      = number of columns of submatrix                            
    !          * if either of nr or nc is nonpositive the code will quit.   
    !                                                                       
    ! ao,                                                                   
    ! jao,iao = extracted matrix in general sparse format with jao containin
    !           the column indices,and iao being the pointer to the beginning  
    !           of the row,in arrays a,ja.                                     
    !----------------------------------------------------------------------c
    !           Y. Saad, Sep. 21 1989                                      c
    !----------------------------------------------------------------------c
    integer i, ii, j, k, k1, k2, klen

    nr = i2-i1+1 
    nc = j2-j1+1 
                                                                       
    if ( nr  <=  0 .or. nc  <=  0) return 
                                                                       
    klen = 0 
                                                                       
    ! Simple procedure. proceeds row-wise...                            
                                                                       
    do 100 i = 1,nr 
      ii = i1+i-1 
      k1 = ia(ii) 
      k2 = ia(ii+1)-1 
      iao(i) = klen+1 
!-----------------------------------------------------------------------
        do 60 k=k1,k2 
          j = ja(k) 
          if (j  >=  j1 .and. j  <=  j2) then 
            klen = klen+1 
            if (job  ==  1) ao(klen) = a(k) 
              jao(klen) = j - j1+1 
            endif 
        60 continue 
    100 continue 
    iao(nr+1) = klen+1 
    return 
    
  end subroutine submat 

!-----------------------------------------------------------------------

      subroutine filter(n,job,drptol,a,ja,ia,b,jb,ib,len,ierr) 
      real(wp) a(*),b(*),drptol 
      integer ja(*),jb(*),ia(*),ib(*),n,job,len,ierr 
!-----------------------------------------------------------------------
!     This module removes any elements whose absolute value             
!     is small from an input matrix A and puts the resulting            
!     matrix in B.  The input parameter job selects a definition        
!     of small.                                                         
!-----------------------------------------------------------------------
! on entry:                                                             
!---------                                                              
!  n         = integer. row dimension of matrix                         
!  job   = integer. used to determine strategy chosen by caller to      
!         drop elements from matrix A.                                  
!          job = 1                                                      
!              Elements whose absolute value is less than the           
!              drop tolerance are removed.                              
!          job = 2                                                      
!              Elements whose absolute value is less than the           
!              product of the drop tolerance and the Euclidean          
!              norm of the row are removed.                             
!          job = 3                                                      
!              Elements whose absolute value is less that the           
!              product of the drop tolerance and the largest            
!              element in the row are removed.                          
!                                                                       
! drptol = real. drop tolerance used for dropping strategy.             
! a                                                                     
! ja                                                                    
! ia     = input matrix in compressed sparse format                     
! len         = integer. the amount of space available in arrays b and j
!                                                                       
! on return:                                                            
!----------                                                             
! b                                                                     
! jb                                                                    
! ib    = resulting matrix in compressed sparse format.                 
!                                                                       
! ierr        = integer. containing error message.                      
!         ierr  ==  0 indicates normal return                           
!         ierr  >  0 indicates that there is'nt enough                 
!         space is a and ja to store the resulting matrix.              
!         ierr then contains the row number where filter stopped.       
! note:                                                                 
!------ This module is in place. (b,jb,ib can ne the same as            
!       a, ja, ia in which case the result will be overwritten).        
!----------------------------------------------------------------------c
!           contributed by David Day,  Sep 19, 1989.                   c
!----------------------------------------------------------------------c
! local variables                                                       
      real(wp) norm,loctol 
      integer index,row,k,k1,k2 
!                                                                       
      index = 1 
      do 10 row= 1,n 
         k1 = ia(row) 
         k2 = ia(row+1) - 1 
         ib(row) = index 
         goto (100,200,300) job 
  100    norm = 1.0d0 
         goto 400 
  200    norm = 0.0d0 
         do 22 k = k1,k2 
            norm = norm + a(k) * a(k) 
   22    continue 
         norm = sqrt(norm) 
         goto 400 
  300    norm = 0.0d0 
         do 23 k = k1,k2 
            if( abs(a(k))   >  norm) then 
               norm = abs(a(k)) 
            endif 
   23    continue 
  400    loctol = drptol * norm 
         do 30 k = k1,k2 
            if( abs(a(k))  >  loctol)then 
               if (index  >  len) then 
               ierr = row 
               return 
            endif 
            b(index) =  a(k) 
            jb(index) = ja(k) 
            index = index + 1 
         endif 
   30 continue 
   10 continue 
      ib(n+1) = index 
      return 
      end subroutine  filter

!-----------------------------------------------------------------------

      subroutine filterm (n,job,drop,a,ja,b,jb,len,ierr) 
      real(wp) a(*),b(*),drop 
      integer ja(*),jb(*),n,job,len,ierr 
!-----------------------------------------------------------------------
!     This subroutine removes any elements whose absolute value         
!     is small from an input matrix A. Same as filter but               
!     uses the MSR format.                                              
!-----------------------------------------------------------------------
! on entry:                                                             
!---------                                                              
!  n         = integer. row dimension of matrix                         
!  job   = integer. used to determine strategy chosen by caller to      
!         drop elements from matrix A.                                  
!          job = 1                                                      
!              Elements whose absolute value is less than the           
!              drop tolerance are removed.                              
!          job = 2                                                      
!              Elements whose absolute value is less than the           
!              product of the drop tolerance and the Euclidean          
!              norm of the row are removed.                             
!          job = 3                                                      
!              Elements whose absolute value is less that the           
!              product of the drop tolerance and the largest            
!              element in the row are removed.                          
!                                                                       
! drop = real. drop tolerance used for dropping strategy.               
! a                                                                     
! ja     = input matrix in Modifief Sparse Row format                   
! len         = integer. the amount of space in arrays b and jb.        
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! b, jb = resulting matrix in Modifief Sparse Row format                
!                                                                       
! ierr        = integer. containing error message.                      
!         ierr  ==  0 indicates normal return                           
!         ierr  >  0 indicates that there is'nt enough                 
!         space is a and ja to store the resulting matrix.              
!         ierr then contains the row number where filter stopped.       
! note:                                                                 
!------ This module is in place. (b,jb can ne the same as               
!       a, ja in which case the result will be overwritten).            
!----------------------------------------------------------------------c
!           contributed by David Day,  Sep 19, 1989.                   c
!----------------------------------------------------------------------c
! local variables                                                       
!                                                                       
      real(wp) norm,loctol 
      integer index,row,k,k1,k2 
!                                                                       
      index = n+2 
      do 10 row= 1,n 
         k1 = ja(row) 
         k2 = ja(row+1) - 1 
         jb(row) = index 
         goto (100,200,300) job 
  100    norm = 1.0d0 
         goto 400 
  200    norm = a(row)**2 
         do 22 k = k1,k2 
            norm = norm + a(k) * a(k) 
   22    continue 
         norm = sqrt(norm) 
         goto 400 
  300    norm = abs(a(row)) 
         do 23 k = k1,k2 
            norm = max(abs(a(k)),norm) 
   23    continue 
  400    loctol = drop * norm 
         do 30 k = k1,k2 
            if( abs(a(k))  >  loctol)then 
               if (index  >  len) then 
                  ierr = row 
                  return 
               endif 
               b(index) =  a(k) 
               jb(index) = ja(k) 
               index = index + 1 
            endif 
   30    continue 
   10 continue 
      jb(n+1) = index 
      return 
      end subroutine  filterm

!-----------------------------------------------------------------------
  subroutine csort ( n, a, ja, ia, iwork, values )

  !*****************************************************************************80
  !
  !! CSORT sorts the elements of a CSR matrix.
  !
  !  Discussion:
  !
  !    This routine sorts the elements of a CSR matrix (stored in Compressed
  !    Sparse Row Format) in increasing order of their column indices within
  !    each row. It uses a form of bucket sort with a cost of O(nnz) where
  !    nnz = number of nonzero elements.
  !
  !    Requires an integer ( kind = 4 ) work array of size length 2*nnz.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! iwork = integer ( kind = 4 ) work array of length max ( n+1, 2*nnz )
  !         where nnz = 2* (ia(n+1)-ia(1))  ) .
  !
  ! values= logical indicating whether or not the real values a(*) must
  !         also be permuted. if (.not. values) then the array a is not
  !         touched by csort and can be a dummy array.
  !
  ! on return:
  !
  ! the matrix stored in the structure a, ja, ia is permuted in such a
  ! way that the column indices are in increasing order within each row.
  ! iwork(1:nnz) contains the permutation used  to rearrange the elements.
  !
    implicit none

    integer ::  n

    real(wp) :: a(*)
    integer :: i
    integer :: ia(n+1)
    integer :: ifirst
    integer :: irow
    integer :: iwork(*)
    integer :: j
    integer :: ja(*)
    integer :: k
    integer :: ko
    integer :: next
    integer :: nnz
    logical :: values
  !
  !  Count the number of elements in each column.
  !
    iwork(1:n+1) = 0

    do i = 1, n
       do k = ia(i), ia(i+1)-1
          j = ja(k) + 1
          iwork(j) = iwork(j) + 1
       end do
    end do
  !
  !  Compute pointers from lengths.
  !
    iwork(1) = 1

    do i = 1, n
       iwork(i+1) = iwork(i) + iwork(i+1)
    end do
  !
  !  Get the positions of the nonzero elements in order of columns.
  !
    ifirst = ia(1)
    nnz = ia(n+1)-ifirst

    do i = 1, n
      do k = ia(i), ia(i+1)-1
        j = ja(k)
        next = iwork(j)
        iwork(nnz+next) = k
        iwork(j) = next + 1
      end do
    end do
  !
  !  Convert to coordinate format.
  !
    do i = 1, n
      do k = ia(i), ia(i+1)-1
        iwork(k) = i
      end do
    end do
  !
  !  Loop to find permutation: for each element find the correct
  !  position in (sorted) arrays A, JA.  Record this in IWORK.
  !
    do k = 1, nnz
       ko = iwork(nnz+k)
       irow = iwork(ko)
       next = ia(irow)
  !
  !  The current element should go in next position in row. IWORK
  !  records this position.
  !
       iwork(ko) = next
       ia(irow) = next + 1
    end do
  !
  !  Perform an in-place permutation of the arrays.
  !
       call ivperm ( nnz, ja(ifirst), iwork )

       if ( values ) then
         call dvperm ( nnz, a(ifirst), iwork )
       end if
  !
  !  Reshift the pointers of the original matrix back.
  !
    do i = n, 1, -1
      ia(i+1) = ia(i)
    end do

    ia(1) = ifirst

    return
  end subroutine csort

!-----------------------------------------------------------------------

      subroutine clean_CSR(nrow, ja,ia) 

      implicit none

!     .. Scalar Arguments ..                                            
      integer, intent(in)   :: nrow
!     ..                                                                
!     .. Array Arguments ..                                             
      integer,  dimension(nrow+1), intent(inout) :: ia
      integer,  dimension(:)     , intent(inout) :: ja

!                                                                       
!     This routine performs two tasks to clean up a CSR matrix          
!     -- remove duplicate entries
!     -- sorts the entries in column increasing order 
!                                                                       
!     nrow    -- row dimension of the matrix                            
!     ncol    -- col dimension of the matrix                            
!     ja,ia   -- input matrix in CSR format                             
!                                                                       
!     On return:                                                        

!     ja,ia   -- compressed and order matrix in CSR format                             

!     local work arrays

      integer                              :: k1
      integer                              :: i,j,k,row_nnz, tot_nnz, del, del_old

      continue

      !if (nrow .le. 0) then
      !  write(*,*) 'Here, nrow', nrow
      !  stop
      !endif

      row_nnz = 0 ; tot_nnz = 0  
      del     = 0 ; del_old = 0

      do i = 1,nrow

        row_nnz = ia(i+1)-ia(i)

        !  Bubble sort the entries

        select case (row_nnz)

          case (0)
            del = row_nnz

          case (1)

            del = row_nnz

            tot_nnz = tot_nnz + 1
            ja(tot_nnz) = ja(ia(i))

          case default

            do k1 = ia(i), ia(i+1)-1
              do j = ia(i), ia(i+1)-2
                if(ja(j+1) < ja(j)) then
                  k       = ja(j+1)
                  ja(j+1) = ja(j)
                  ja(j)   = k 
                endif
              enddo
            enddo

            del = 1

            tot_nnz = tot_nnz + 1
            ja(tot_nnz) = ja(ia(i))

            do j = ia(i)+1, ia(i+1)-1
              if(ja(j) /= ja(j-1)) then
                tot_nnz = tot_nnz + 1
                del   = del   + 1
                ja(tot_nnz) = ja(j) 
              endif
            enddo

        end select

        if(i /= 1) then 
          ia(i) = ia(i-1) + del_old
          del_old = del
        else
          del_old = del
        endif

      enddo

      ia(nrow+1) = ia(nrow) + del_old

      end subroutine clean_CSR

!-----------------------------------------------------------------------

      subroutine clncsr(job,value2,nrow,a,ja,ia,indu,iwk) 
!     .. Scalar Arguments ..                                            
      integer job, nrow, value2 
!     ..                                                                
!     .. Array Arguments ..                                             
      integer ia(nrow+1),indu(nrow),iwk(nrow+1),ja(*) 
      real(wp)  a(*) 
!     ..                                                                
!                                                                       
!     This routine performs two tasks to clean up a CSR matrix          
!     -- remove duplicate/zero entries,                                 
!     -- perform a partial ordering, new order lower triangular part,   
!        main diagonal, upper triangular part.                          
!                                                                       
!     On entry:                                                         
!                                                                       
!     job   = options                                                   
!         0 -- nothing is done                                          
!         1 -- eliminate duplicate entries, zero entries.               
!         2 -- eliminate duplicate entries and perform partial ordering.
!         3 -- eliminate duplicate entries, sort the entries in the     
!              increasing order of column indices.                       
!                                                                       
!     value2  -- 0 the matrix is pattern only (a is not touched)        
!                1 matrix has values too.                               
!     nrow    -- row dimension of the matrix                            
!     a,ja,ia -- input matrix in CSR format                             
!                                                                       
!     On return:                                                        
!     a,ja,ia -- cleaned matrix                                         
!     indu    -- pointers to the beginning of the upper triangular      
!                portion if job > 1                                     
!                                                                       
!     Work space:                                                       
!     iwk     -- integer work space of size nrow+1                      
!                                                                       
!     .. Local Scalars ..                                               
      integer i,j,k,ko,ipos,kfirst,klast 
      real(wp)  tmp 
!     ..                                                                
!                                                                       
      if (job <= 0) return 
!                                                                       
!     .. eliminate duplicate entries --                                 
!     array INDU is used as marker for existing indices, it is also the 
!     location of the entry.                                            
!     IWK is used to stored the old IA array.                           
!     matrix is copied to squeeze out the space taken by the duplicated 
!     entries.                                                          
!                                                                       
      do 90 i = 1, nrow 
         indu(i) = 0 
         iwk(i) = ia(i) 
   90 continue 
      iwk(nrow+1) = ia(nrow+1) 
      k = 1 
      do 120 i = 1, nrow 
         ia(i) = k 
         ipos = iwk(i) 
         klast = iwk(i+1) 
  100    if (ipos < klast) then 
            j = ja(ipos) 
            if (indu(j) == 0) then 
!     .. new entry ..                                                   
               if (value2 /= 0) then 
                  if (a(ipos)  /=  0.0D0) then 
                     indu(j) = k 
                     ja(k) = ja(ipos) 
                     a(k) = a(ipos) 
                     k = k + 1 
                  endif 
               else 
                  indu(j) = k 
                  ja(k) = ja(ipos) 
                  k = k + 1 
               endif 
            else if (value2 /= 0) then 
!     .. duplicate entry ..                                             
               a(indu(j)) = a(indu(j)) + a(ipos) 
            endif 
            ipos = ipos + 1 
            go to 100 
         endif 
!     .. remove marks before working on the next row ..                 
         do 110 ipos = ia(i), k - 1 
            indu(ja(ipos)) = 0 
  110    continue 
  120 continue 
      ia(nrow+1) = k 
      if (job <= 1) return 
!                                                                       
!     .. partial ordering ..                                            
!     split the matrix into strict upper/lower triangular               
!     parts, INDU points to the the beginning of the upper part.        
!                                                                       
      do 140 i = 1, nrow 
         klast = ia(i+1) - 1 
         kfirst = ia(i) 
  130    if (klast > kfirst) then 
            if (ja(klast) < i .and. ja(kfirst) >= i) then 
!     .. swap klast with kfirst ..                                      
               j = ja(klast) 
               ja(klast) = ja(kfirst) 
               ja(kfirst) = j 
               if (value2 /= 0) then 
                  tmp = a(klast) 
                  a(klast) = a(kfirst) 
                  a(kfirst) = tmp 
               endif 
            endif 
            if (ja(klast) >= i)                                         &
     &         klast = klast - 1                                        
            if (ja(kfirst) < i)                                        &
     &         kfirst = kfirst + 1                                      
            go to 130 
         endif 
!                                                                       
         if (ja(klast) < i) then 
            indu(i) = klast + 1 
         else 
            indu(i) = klast 
         endif 
  140 continue 
      if (job <= 2) return 
!                                                                       
!     .. order the entries according to column indices                  
!     burble-sort is used                                               
!                                                                       
      do 190 i = 1, nrow 
         do 160 ipos = ia(i), indu(i)-1 
            do 150 j = indu(i)-1, ipos+1, -1 
               k = j - 1 
               if (ja(k) > ja(j)) then 
                  ko = ja(k) 
                  ja(k) = ja(j) 
                  ja(j) = ko 
                  if (value2 /= 0) then 
                     tmp = a(k) 
                     a(k) = a(j) 
                     a(j) = tmp 
                  endif 
               endif 
  150       continue 
  160    continue 
         do 180 ipos = indu(i), ia(i+1)-1 
            do 170 j = ia(i+1)-1, ipos+1, -1 
               k = j - 1 
               if (ja(k) > ja(j)) then 
                  ko = ja(k) 
                  ja(k) = ja(j) 
                  ja(j) = ko 
                  if (value2 /= 0) then 
                     tmp = a(k) 
                     a(k) = a(j) 
                     a(j) = tmp 
                  endif 
               endif 
  170       continue 
  180    continue 
  190 continue 
      return 
      end subroutine  clncsr

!-----------------------------------------------------------------------

      subroutine copmat (nrow,a,ja,ia,ao,jao,iao,ipos,job) 
      real(wp) a(*),ao(*) 
      integer nrow, ia(*),ja(*),jao(*),iao(*), ipos, job 
!---------------------------------------------------------------------- 
! copies the matrix a, ja, ia, into the matrix ao, jao, iao.            
!---------------------------------------------------------------------- 
! on entry:                                                             
!---------                                                              
! nrow        = row dimension of the matrix                             
! a,                                                                    
! ja,                                                                   
! ia    = input matrix in compressed sparse row format.                 
! ipos  = integer. indicates the position in the array ao, jao          
!         where the first element should be copied. Thus                
!         iao(1) = ipos on return.                                      
! job   = job indicator. if (job  /=  1) the values are not copies      
!         (i.e., pattern only is copied in the form of arrays ja, ia).  
!                                                                       
! on return:                                                            
!----------                                                             
! ao,                                                                   
! jao,                                                                  
! iao   = output matrix containing the same data as a, ja, ia.          
!-----------------------------------------------------------------------
!           Y. Saad, March 1990.                                        
!-----------------------------------------------------------------------
! local variables                                                       
      integer kst, i, k 
!                                                                       
      kst    = ipos -ia(1) 
      do 100 i = 1, nrow+1 
         iao(i) = ia(i) + kst 
  100 continue 
!                                                                       
      do 200 k=ia(1), ia(nrow+1)-1 
         jao(kst+k)= ja(k) 
  200 continue 
!                                                                       
      if (job  /=  1) return 
      do 201 k=ia(1), ia(nrow+1)-1 
         ao(kst+k) = a(k) 
  201 continue 
!                                                                       
      return 
      end subroutine  copmat

!-----------------------------------------------------------------------

      subroutine msrcop (nrow,a,ja,ao,jao,job) 
      real(wp) a(*),ao(*) 
      integer nrow, ja(*),jao(*), job 
!---------------------------------------------------------------------- 
! copies the MSR matrix a, ja, into the MSR matrix ao, jao              
!---------------------------------------------------------------------- 
! on entry:                                                             
!---------                                                              
! nrow        = row dimension of the matrix                             
! a,ja  = input matrix in Modified compressed sparse row format.        
! job   = job indicator. Values are not copied if job  /=  1            
!                                                                       
! on return:                                                            
!----------                                                             
! ao, jao   = output matrix containing the same data as a, ja.          
!-----------------------------------------------------------------------
!           Y. Saad,                                                    
!-----------------------------------------------------------------------
! local variables                                                       
      integer i, k 
!                                                                       
      do 100 i = 1, nrow+1 
         jao(i) = ja(i) 
  100 continue 
!                                                                       
      do 200 k=ja(1), ja(nrow+1)-1 
         jao(k)= ja(k) 
  200 continue 
!                                                                       
      if (job  /=  1) return 
      do 201 k=ja(1), ja(nrow+1)-1 
         ao(k) = a(k) 
  201 continue 
      do 202 k=1,nrow 
         ao(k) = a(k) 
  202 continue 
!                                                                       
      return 
      end subroutine  msrcop

!-----------------------------------------------------------------------

      double precision function getelm (i,j,a,ja,ia,iadd,sorted) 
!-----------------------------------------------------------------------
!     purpose:                                                          
!     --------                                                          
!     this function returns the element a(i,j) of a matrix a,           
!     for any pair (i,j).  the matrix is assumed to be stored           
!     in compressed sparse row (csr) format. getelm performs a          
!     binary search in the case where it is known that the elements     
!     are sorted so that the column indices are in increasing order.    
!     also returns (in iadd) the address of the element a(i,j) in       
!     arrays a and ja when the search is successsful (zero if not).     
!-----                                                                  
!     first contributed by noel nachtigal (mit).                        
!     recoded jan. 20, 1991, by y. saad [in particular                  
!     added handling of the non-sorted case + the iadd output]          
!-----------------------------------------------------------------------
!     parameters:                                                       
!     -----------                                                       
! on entry:                                                             
!----------                                                             
!     i      = the row index of the element sought (input).             
!     j      = the column index of the element sought (input).          
!     a      = the matrix a in compressed sparse row format (input).    
!     ja     = the array of column indices (input).                     
!     ia     = the array of pointers to the rows' data (input).         
!     sorted = logical indicating whether the matrix is knonw to        
!              have its column indices sorted in increasing order       
!              (sorted=.true.) or not (sorted=.false.).                 
!              (input).                                                 
! on return:                                                            
!-----------                                                            
!     getelm = value of a(i,j).                                         
!     iadd   = address of element a(i,j) in arrays a, ja if found,      
!              zero if not found. (output)                              
!                                                                       
!     note: the inputs i and j are not checked for validity.            
!-----------------------------------------------------------------------
!     noel m. nachtigal october 28, 1990 -- youcef saad jan 20, 1991.   
!-----------------------------------------------------------------------
      integer i, ia(*), iadd, j, ja(*) 
      double precision a(*) 
      logical sorted 
!                                                                       
!     local variables.                                                  
!                                                                       
      integer ibeg, iend, imid, k 
!                                                                       
!     initialization                                                    
!                                                                       
      iadd = 0 
      getelm = 0.0 
      ibeg = ia(i) 
      iend = ia(i+1)-1 
!                                                                       
!     case where matrix is not necessarily sorted                       
!                                                                       
      if (.not. sorted) then 
!                                                                       
! scan the row - exit as soon as a(i,j) is found                        
!                                                                       
         do 5  k=ibeg, iend 
            if (ja(k)  ==   j) then 
               iadd = k 
               goto 20 
            endif 
    5    continue 
!                                                                       
!     end unsorted case. begin sorted case                              
!                                                                       
      else 
!                                                                       
!     begin binary search.   compute the middle index.                  
!                                                                       
   10    imid = ( ibeg + iend ) / 2 
!                                                                       
!     test if  found                                                    
!                                                                       
         if (ja(imid) == j) then 
            iadd = imid 
            goto 20 
         endif 
         if (ibeg  >=  iend) goto 20 
!                                                                       
!     else     update the interval bounds.                              
!                                                                       
         if (ja(imid) > j) then 
            iend = imid -1 
         else 
            ibeg = imid +1 
         endif 
         goto 10 
!                                                                       
!     end both cases                                                    
!                                                                       
      endif 
!                                                                       
   20 if (iadd  /=  0) getelm = a(iadd) 
!                                                                       
      return 
      end function getelm

!-----------------------------------------------------------------------

      subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff) 
      real(wp) diag(*),a(*) 
      integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*) 
!-----------------------------------------------------------------------
! this subroutine extracts a given diagonal from a matrix stored in csr 
! format. the output matrix may be transformed with the diagonal removed
! from it if desired (as indicated by job.)                             
!-----------------------------------------------------------------------
! our definition of a diagonal of matrix is a vector of length nrow     
! (always) which contains the elements in rows 1 to nrow of             
! the matrix that are contained in the diagonal offset by ioff          
! with respect to the main diagonal. if the diagonal element            
! falls outside the matrix then it is defined as a zero entry.          
! thus the proper definition of diag(*) with offset ioff is             
!                                                                       
!     diag(i) = a(i,ioff+i) i=1,2,...,nrow                              
!     with elements falling outside the matrix being defined as zero.   
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
! on entry:                                                             
!----------                                                             
!                                                                       
! nrow        = integer. the row dimension of the matrix a.             
! ncol        = integer. the column dimension of the matrix a.          
! job   = integer. job indicator.  if job = 0 then                      
!         the matrix a, ja, ia, is not altered on return.               
!         if job /= 0  then getdia will remove the entries              
!         collected in diag from the original matrix.                   
!         this is done in place.                                        
!                                                                       
! a,ja,                                                                 
!    ia = matrix stored in compressed sparse row a,ja,ia,format         
! ioff  = integer,containing the offset of the wanted diagonal          
!          the diagonal extracted is the one corresponding to the       
!          entries a(i,j) with j-i = ioff.                              
!          thus ioff = 0 means the main diagonal                        
!                                                                       
! on return:                                                            
!-----------                                                            
! len   = number of nonzero elements found in diag.                     
!         (len  <=  min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )             
!                                                                       
! diag  = real(wp) array of length nrow containing the wanted diagonal.   
!          diag contains the diagonal (a(i,j),j-i = ioff ) as defined   
!         above.                                                        
!                                                                       
! idiag = integer array of  length len, containing the poisitions       
!         in the original arrays a and ja of the diagonal elements      
!         collected in diag. a zero entry in idiag(i) means that        
!         there was no entry found in row i belonging to the diagonal.  
!                                                                       
! a, ja,                                                                
!    ia = if job  /=  0 the matrix is unchanged. otherwise the nonzero  
!         diagonal entries collected in diag are removed from the       
!         matrix and therefore the arrays a, ja, ia will change.        
!          (the matrix a, ja, ia will contain len fewer elements)       
!                                                                       
!----------------------------------------------------------------------c
!     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c
!----------------------------------------------------------------------c
!     local variables                                                   
      integer istart, max, iend, i, kold, k, kdiag, ko 
!                                                                       
      istart = max(0,-ioff) 
      iend = min(nrow,ncol-ioff) 
      len = 0 
      do 1 i=1,nrow 
         idiag(i) = 0 
         diag(i) = 0.0d0 
    1 continue 
!                                                                       
!     extract  diagonal elements                                        
!                                                                       
      do 6 i=istart+1, iend 
         do 51 k= ia(i),ia(i+1) -1 
            if (ja(k)-i  ==  ioff) then 
               diag(i)= a(k) 
               idiag(i) = k 
               len = len+1 
               goto 6 
            endif 
   51    continue 
    6 continue 
      if (job  ==  0 .or. len  == 0) return 
!                                                                       
!     remove diagonal elements and rewind structure                     
!                                                                       
      ko = 0 
      do  7 i=1, nrow 
         kold = ko 
         kdiag = idiag(i) 
         do 71 k= ia(i), ia(i+1)-1 
            if (k  /=  kdiag) then 
               ko = ko+1 
               a(ko) = a(k) 
               ja(ko) = ja(k) 
            endif 
   71    continue 
         ia(i) = kold+1 
    7 continue 
!                                                                       
!     redefine ia(nrow+1)                                               
!                                                                       
      ia(nrow+1) = ko+1 
      return 
      end subroutine  getdia

!-----------------------------------------------------------------------

      subroutine transp (nrow,ncol,a,ja,ia,iwk,ierr) 
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr 
      real(wp) a(*) 
!-----------------------------------------------------------------------
! In-place transposition routine.                                       
!-----------------------------------------------------------------------
! this subroutine transposes a matrix stored in compressed sparse row   
! format. the transposition is done in place in that the arrays a,ja,ia 
! of the transpose are overwritten onto the original arrays.            
!-----------------------------------------------------------------------
! on entry:                                                             
!---------                                                              
! nrow        = integer. The row dimension of A.                        
! ncol        = integer. The column dimension of A.                     
! a        = real array of size nnz (number of nonzero elements in A).  
!         containing the nonzero elements                               
! ja        = integer array of length nnz containing the column position
!           of the corresponding elements in a.                         
! ia        = integer of size n+1, where n = max(nrow,ncol). On entry   
!         ia(k) contains the position in a,ja of  the beginning of      
!         the k-th row.                                                 
!                                                                       
! iwk        = integer work array of same length as ja.                 
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! ncol        = actual row dimension of the transpose of the input matri
!         Note that this may be  <=  the input value for ncol, in       
!         case some of the last columns of the input matrix are zero    
!         columns. In the case where the actual number of rows found    
!         in transp(A) exceeds the input value of ncol, transp will     
!         return without completing the transposition. see ierr.        
! a,                                                                    
! ja,                                                                   
! ia        = contains the transposed matrix in compressed sparse       
!         row format. The row dimension of a, ja, ia is now ncol.       
!                                                                       
! ierr        = integer. error message. If the number of rows for the   
!         transposed matrix exceeds the input value of ncol,            
!         then ierr is  set to that number and transp quits.            
!         Otherwise ierr is set to 0 (normal return).                   
!                                                                       
! Note:                                                                 
!----- 1) If you do not need the transposition to be done in place      
!         it is preferrable to use the conversion routine csrcsc        
!         (see conversion routines in formats).                         
!      2) the entries of the output matrix are not sorted (the column   
!         indices in each are not in increasing order) use csrcsc       
!         if you want them sorted.                                      
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!  modified Oct. 11, 1989.                                             c
!----------------------------------------------------------------------c
! local variables                                                       
      real(wp) t, t1 
      integer i, inext, init, j, jcol, k, nnz, l
      ierr = 0 
      nnz = ia(nrow+1)-1 
!                                                                       
!     determine column dimension                                        
!                                                                       
      jcol = 0 
      do 1 k=1, nnz 
         jcol = max(jcol,ja(k)) 
    1 continue 
      if (jcol  >  ncol) then 
         ierr = jcol 
         return 
      endif 
!                                                                       
!     convert to coordinate format. use iwk for row indices.            
!                                                                       
      ncol = jcol 
!                                                                       
      do 3 i=1,nrow 
         do 2 k=ia(i),ia(i+1)-1 
            iwk(k) = i 
    2    continue 
    3 continue 
!     find pointer array for transpose.                                 
      do 35 i=1,ncol+1 
         ia(i) = 0 
   35 continue 
      do 4 k=1,nnz 
         i = ja(k) 
         ia(i+1) = ia(i+1)+1 
    4 continue 
      ia(1) = 1 
!-----------------------------------------------------------------------
      do 44 i=1,ncol 
         ia(i+1) = ia(i) + ia(i+1) 
   44 continue 
!                                                                       
!     loop for a cycle in chasing process.                              
!                                                                       
      init = 1 
      k = 0 
    5 t = a(init) 
      i = ja(init) 
      j = iwk(init) 
      iwk(init) = -1 
!-----------------------------------------------------------------------
    6 k = k+1 
!     current row number is i.  determine  where to go.                 
      l = ia(i) 
!     save the chased element.                                          
      t1 = a(l) 
      inext = ja(l) 
!     then occupy its location.                                         
      a(l)  = t 
      ja(l) = j 
!     update pointer information for next element to be put in row i.   
      ia(i) = l+1 
!     determine  next element to be chased                              
      if (iwk(l)  <  0) goto 65 
      t = t1 
      i = inext 
      j = iwk(l) 
      iwk(l) = -1 
      if (k  <  nnz) goto 6 
      goto 70 
   65 init = init+1 
      if (init  >  nnz) goto 70 
      if (iwk(init)  <  0) goto 65 
!     restart chasing --                                                
      goto 5 
   70 continue 
      do 80 i=ncol,1,-1 
         ia(i+1) = ia(i) 
   80 continue 
      ia(1) = 1 
!                                                                       
      return 
      end subroutine  transp

!-----------------------------------------------------------------------

      subroutine getl (n,a,ja,ia,ao,jao,iao) 
      integer n, ia(*), ja(*), iao(*), jao(*) 
      real(wp) a(*), ao(*) 
!-----------------------------------------------------------------------
! this subroutine extracts the lower triangular part of a matrix        
! and writes the result ao, jao, iao. The routine is in place in        
! that ao, jao, iao can be the same as a, ja, ia if desired.            
!-----------                                                            
! on input:                                                             
!                                                                       
! n     = dimension of the matrix a.                                    
! a, ja,                                                                
!    ia = matrix stored in compressed sparse row format.                
! On return:                                                            
! ao, jao,                                                              
!    iao = lower triangular matrix (lower part of a)                    
!        stored in a, ja, ia, format                                    
! note: the diagonal element is the last element in each row.           
! i.e. in  a(ia(i+1)-1 )                                                
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case   
! getl will overwrite the result on a, ja, ia.                          
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
      real(wp) t 
      integer ko, kold, kdiag, k, i 
!                                                                       
! inititialize ko (pointer for output matrix)                           
!                                                                       
      ko = 0 
      do  7 i=1, n 
         kold = ko 
         kdiag = 0 
         do 71 k = ia(i), ia(i+1) -1 
            if (ja(k)   >  i) goto 71 
            ko = ko+1 
            ao(ko) = a(k) 
            jao(ko) = ja(k) 
            if (ja(k)   ==  i) kdiag = ko 
   71    continue 
         if (kdiag  ==  0 .or. kdiag  ==  ko) goto 72 
!                                                                       
!     exchange                                                          
!                                                                       
         t = ao(kdiag) 
         ao(kdiag) = ao(ko) 
         ao(ko) = t 
!                                                                       
         k = jao(kdiag) 
         jao(kdiag) = jao(ko) 
         jao(ko) = k 
   72    iao(i) = kold+1 
    7 continue 
!     redefine iao(n+1)                                                 
      iao(n+1) = ko+1 
      return 
      end subroutine  getl

!-----------------------------------------------------------------------

      subroutine getu (n,a,ja,ia,ao,jao,iao) 
      integer n, ia(*), ja(*), iao(*), jao(*) 
      real(wp) a(*), ao(*) 
!-----------------------------------------------------------------------
! this subroutine extracts the upper triangular part of a matrix        
! and writes the result ao, jao, iao. The routine is in place in        
! that ao, jao, iao can be the same as a, ja, ia if desired.            
!-----------                                                            
! on input:                                                             
!                                                                       
! n     = dimension of the matrix a.                                    
! a, ja,                                                                
!    ia = matrix stored in a, ja, ia, format                            
! On return:                                                            
! ao, jao,                                                              
!    iao = upper triangular matrix (upper part of a)                    
!        stored in compressed sparse row format                         
! note: the diagonal element is the last element in each row.           
! i.e. in  a(ia(i+1)-1 )                                                
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case   
! getu will overwrite the result on a, ja, ia.                          
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
      real(wp) t 
      integer ko, k, i, kdiag, kfirst 
      ko = 0 
      do  7 i=1, n 
         kfirst = ko+1 
         kdiag = 0 
         do 71 k = ia(i), ia(i+1) -1 
            if (ja(k)   <  i) goto 71 
            ko = ko+1 
            ao(ko) = a(k) 
            jao(ko) = ja(k) 
            if (ja(k)   ==  i) kdiag = ko 
   71    continue 
         if (kdiag  ==  0 .or. kdiag  ==  kfirst) goto 72 
!     exchange                                                          
         t = ao(kdiag) 
         ao(kdiag) = ao(kfirst) 
         ao(kfirst) = t 
!                                                                       
         k = jao(kdiag) 
         jao(kdiag) = jao(kfirst) 
         jao(kfirst) = k 
   72    iao(i) = kfirst 
    7 continue 
!     redefine iao(n+1)                                                 
      iao(n+1) = ko+1 
      return 
      end subroutine  getu

!-----------------------------------------------------------------------

      subroutine levels (n, jal, ial, nlev, lev, ilev, levnum) 
      integer n, nlev
      integer jal(*),ial(*), levnum(*), ilev(*), lev(*) 
!-----------------------------------------------------------------------
! levels gets the level structure of a lower triangular matrix          
! for level scheduling in the parallel solution of triangular systems   
! strict lower matrices (e.g. unit) as well matrices with their main    
! diagonal are accepted.                                                
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n        = integer. The row dimension of the matrix                   
! jal, ial =                                                            
!                                                                       
! on return:                                                            
!-----------                                                            
! nlev     = integer. number of levels found                            
! lev      = integer array of length n containing the level             
!            scheduling permutation.                                    
! ilev     = integer array. pointer to beginning of levels in lev.      
!            the numbers lev(i) to lev(i+1)-1 contain the row numbers   
!            that belong to level number i, in the level scheduling     
!            ordering. The equations of the same level can be solved    
!            in parallel, once those of all the previous levels have    
!            been solved.                                               
! work arrays:                                                          
!-------------                                                          
! levnum   = integer array of length n (containing the level numbers    
!            of each unknown on return)                                 
!-----------------------------------------------------------------------
      integer i, j, levi    

      do 10 i = 1, n 
         levnum(i) = 0 
   10 continue 
!                                                                       
!     compute level of each node --                                     
!                                                                       
      nlev = 0 
      do 20 i = 1, n 
         levi = 0 
         do 15 j = ial(i), ial(i+1) - 1 
            levi = max (levi, levnum(jal(j))) 
   15    continue 
         levi = levi+1 
         levnum(i) = levi 
         nlev = max(nlev,levi) 
   20 continue 
!-------------set data structure  --------------------------------------
      do 21 j=1, nlev+1 
         ilev(j) = 0 
   21 continue 
!------count  number   of elements in each level -----------------------
      do 22 j=1, n 
         i = levnum(j)+1 
         ilev(i) = ilev(i)+1 
   22 continue 
!---- set up pointer for  each  level ----------------------------------
      ilev(1) = 1 
      do 23 j=1, nlev 
         ilev(j+1) = ilev(j)+ilev(j+1) 
   23 continue 
!-----determine elements of each level -------------------------------- 
      do 30 j=1,n 
         i = levnum(j) 
         lev(ilev(i)) = j 
         ilev(i) = ilev(i)+1 
   30 continue 
!     reset pointers backwards                                          
      do 35 j=nlev, 1, -1 
         ilev(j+1) = ilev(j) 
   35 continue 
      ilev(1) = 1 
      return 
      end subroutine  levels

!-----------------------------------------------------------------------

      subroutine amask (nrow,ncol,a,ja,ia,jmask,imask,                  &
     &                  c,jc,ic,iw,nzmax,ierr)                          
!---------------------------------------------------------------------  
      real(wp) a(*),c(*)
      integer nrow, ncol, nzmax, ierr
      integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*),imask(nrow+1) 
      logical iw(ncol) 
!-----------------------------------------------------------------------
! This subroutine builds a sparse matrix from an input matrix by        
! extracting only elements in positions defined by the mask jmask, imask
!-----------------------------------------------------------------------
! On entry:                                                             
!---------                                                              
! nrow  = integer. row dimension of input matrix                        
! ncol        = integer. Column dimension of input matrix.              
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia        = matrix in Compressed Sparse Row format                    
!                                                                       
! jmask,                                                                
! imask = matrix defining mask (pattern only) stored in compressed      
!         sparse row format.                                            
!                                                                       
! nzmax = length of arrays c and jc. see ierr.                          
!                                                                       
! On return:                                                            
!-----------                                                            
!                                                                       
! a, ja, ia and jmask, imask are unchanged.                             
!                                                                       
! c                                                                     
! jc,                                                                   
! ic        = the output matrix in Compressed Sparse Row format.        
!                                                                       
! ierr  = integer. serving as error message.c                           
!         ierr = 1  means normal return                                 
!         ierr  >  1 means that amask stopped when processing          
!         row number ierr, because there was not enough space in        
!         c, jc according to the value of nzmax.                        
!                                                                       
! work arrays:                                                          
!-------------                                                          
! iw        = logical work array of length ncol.                        
!                                                                       
! note:                                                                 
!------ the  algorithm is in place: c, jc, ic can be the same as        
! a, ja, ia in which cas the code will overwrite the matrix c           
! on a, ja, ia                                                          
!                                                                       
!-----------------------------------------------------------------------
      integer ii, j, k, k1, k2, len

      ierr = 0 
      len = 0 
      do 1 j=1, ncol 
         iw(j) = .false. 
    1 continue 
!     unpack the mask for row ii in iw                                  
      do 100 ii=1, nrow 
!     save pointer in order to be able to do things in place            
         do 2 k=imask(ii), imask(ii+1)-1 
            iw(jmask(k)) = .true. 
    2    continue 
!     add umasked elemnts of row ii                                     
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         ic(ii) = len+1 
         do 200 k=k1,k2 
            j = ja(k) 
            if (iw(j)) then 
               len = len+1 
               if (len  >  nzmax) then 
                  ierr = ii 
                  return 
               endif 
               jc(len) = j 
               c(len) = a(k) 
            endif 
  200    continue 
!                                                                       
         do 3 k=imask(ii), imask(ii+1)-1 
            iw(jmask(k)) = .false. 
    3    continue 
  100 continue 
      ic(nrow+1)=len+1 
!                                                                       
      return 
      end subroutine  amask

!-----------------------------------------------------------------------

      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job 
      real(wp) a(*),ao(*) 
!-----------------------------------------------------------------------
! this subroutine permutes the rows of a matrix in CSR format.          
! rperm  computes B = P A  where P is a permutation matrix.             
! the permutation P is defined through the array perm: for each j,      
! perm(j) represents the destination row number of row number j.        
! Youcef Saad -- recoded Jan 28, 1991.                                  
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n         = dimension of the matrix                                   
! a, ja, ia = input matrix in csr format                                
! perm         = integer array of length nrow containing the permutation
!          for the rows: perm(i) is the destination of row i in the     
!         permuted matrix.                                              
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)       
!         in the output  matrix.                                        
!                                                                       
! job        = integer indicating the work to be done:                  
!                 job = 1        permute a, ja, ia into ao, jao, iao    
!                       (including the copying of real values ao and    
!                       the array iao).                                 
!                 job  /=  1 :  ignore real values.                     
!                     (in which case arrays a and ao are not needed nor 
!                      used).                                           
!                                                                       
!------------                                                           
! on return:                                                            
!------------                                                           
! ao, jao, iao = input matrix in a, ja, ia format                       
! note :                                                                
!        if (job /= 1)  then the arrays a and ao are not used.          
!----------------------------------------------------------------------c
!           Y. Saad, May  2, 1990                                      c
!----------------------------------------------------------------------c
      integer i , ii, j, k, ko
      logical values 
      
      values = (job  ==  1) 
!                                                                       
!     determine pointers for output matix.                              
!                                                                       
      do 50 j=1,nrow 
         i = perm(j) 
         iao(i+1) = ia(j+1) - ia(j) 
   50 continue 
!                                                                       
! get pointers from lengths                                             
!                                                                       
      iao(1) = 1 
      do 51 j=1,nrow 
         iao(j+1)=iao(j+1)+iao(j) 
   51 continue 
!                                                                       
! copying                                                               
!                                                                       
      do 100 ii=1,nrow 
!                                                                       
! old row = ii  -- new row = iperm(ii) -- ko = new pointer              
!                                                                       
         ko = iao(perm(ii)) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k) 
            ko = ko+1 
   60    continue 
  100 continue 
!                                                                       
      return 
      end subroutine  rperm

!-----------------------------------------------------------------------

      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job 
      real(wp) a(*), ao(*) 
!-----------------------------------------------------------------------
! this subroutine permutes the columns of a matrix a, ja, ia.           
! the result is written in the output matrix  ao, jao, iao.             
! cperm computes B = A P, where  P is a permutation matrix              
! that maps column j into column perm(j), i.e., on return               
!      a(i,j) becomes a(i,perm(j)) in new matrix                        
! Y. Saad, May 2, 1990 / modified Jan. 28, 1991.                        
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! nrow         = row dimension of the matrix                            
!                                                                       
! a, ja, ia = input matrix in csr format.                               
!                                                                       
! perm        = integer array of length ncol (number of columns of A    
!         containing the permutation array  the columns:                
!         a(i,j) in the original matrix becomes a(i,perm(j))            
!         in the output matrix.                                         
!                                                                       
! job        = integer indicating the work to be done:                  
!                 job = 1        permute a, ja, ia into ao, jao, iao    
!                       (including the copying of real values ao and    
!                       the array iao).                                 
!                 job  /=  1 :  ignore real values ao and ignore iao.   
!                                                                       
!------------                                                           
! on return:                                                            
!------------                                                           
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed) 
!                                                                       
! Notes:                                                                
!-------                                                                
! 1. if job=1 then ao, iao are not used.                                
! 2. This routine is in place: ja, jao can be the same.                 
! 3. If the matrix is initially sorted (by increasing column number)    
!    then ao,jao,iao  may not be on return.                             
!                                                                       
!----------------------------------------------------------------------c
! local parameters:                                                     
      integer k, i, nnz 
!                                                                       
      nnz = ia(nrow+1)-1 
      do 100 k=1,nnz 
         jao(k) = perm(ja(k)) 
  100 continue 
!                                                                       
!     done with ja array. return if no need to touch values.            
!                                                                       
      if (job  /=  1) return 
!                                                                       
! else get new pointers -- and copy values too.                         
!                                                                       
      do 1 i=1, nrow+1 
         iao(i) = ia(i) 
    1 continue 
!                                                                       
      do 2 k=1, nnz 
         ao(k) = a(k) 
    2 continue 
!                                                                       
      return 
      end subroutine  cperm

!-----------------------------------------------------------------------

      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),      &
     &        qperm(*),job                                              
      real(wp) a(*),ao(*) 
!-----------------------------------------------------------------------
! This routine permutes the rows and columns of a matrix stored in CSR  
! format. i.e., it computes P A Q, where P, Q are permutation matrices. 
! P maps row i into row perm(i) and Q maps column j into column qperm(j)
!      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix            
! In the particular case where Q is the transpose of P (symmetric       
! permutation of A) then qperm is not needed.                           
! note that qperm should be of length ncol (number of columns) but this 
! is not checked.                                                       
!-----------------------------------------------------------------------
! Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.                         
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n         = dimension of the matrix                                   
! a, ja,                                                                
!    ia = input matrix in a, ja, ia format                              
! perm         = integer array of length n containing the permutation ar
!          for the rows: perm(i) is the destination of row i in the     
!         permuted matrix -- also the destination of column i in case   
!         permutation is symmetric (job  <=  2)                         
!                                                                       
! qperm        = same thing for the columns. This should be provided onl
!         if job=3 or job=4, i.e., only in the case of a nonsymmetric   
!          permutation of rows and columns. Otherwise qperm is a dummy  
!                                                                       
! job        = integer indicating the work to be done:                  
! * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)        
!                 job = 1        permute a, ja, ia into ao, jao, iao    
!                 job = 2 permute matrix ignoring real values.          
! * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q            
!                 job = 3        permute a, ja, ia into ao, jao, iao    
!                 job = 4 permute matrix ignoring real values.          
!                                                                       
! on return:                                                            
!-----------                                                            
! ao, jao, iao = input matrix in a, ja, ia format                       
!                                                                       
! in case job  ==  2 or job  ==  4, a and ao are never referred to      
! and can be dummy arguments.                                           
! Notes:                                                                
!-------                                                                
!  1) algorithm is in place                                             
!  2) column indices may not be sorted on return even  though they may b
!     on entry.                                                         
!----------------------------------------------------------------------c
! local variables                                                       
      integer locjob, mod 
!                                                                       
!     locjob indicates whether or not real values must be copied.       
!                                                                       
      locjob = mod(job,2) 
!                                                                       
! permute rows first                                                    
!                                                                       
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob) 
!                                                                       
! then permute columns                                                  
!                                                                       
      locjob = 0 
!                                                                       
      if (job  <=  2) then 
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob) 
      else 
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob) 
      endif 
!                                                                       
      return 
      end subroutine  dperm

!-----------------------------------------------------------------------

      subroutine dperm1 (i1,i2,a,ja,ia,b,jb,ib,perm,ipos,job) 
      integer ipos
      integer i1,i2,job,ja(*),ia(*),jb(*),ib(*),perm(*) 
      real(wp) a(*),b(*) 
!-----------------------------------------------------------------------
!     general submatrix extraction routine.                             
!-----------------------------------------------------------------------
!     extracts rows perm(i1), perm(i1+1), ..., perm(i2) (in this order) 
!     from a matrix (doing nothing in the column indices.) The resulting
!     submatrix is constructed in b, jb, ib. A pointer ipos to the      
!     beginning of arrays b,jb,is also allowed (i.e., nonzero elements  
!     are accumulated starting in position ipos of b, jb).              
!-----------------------------------------------------------------------
! Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991 / modified for PSPARSLIB  
! Sept. 1997..                                                          
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n         = dimension of the matrix                                   
! a,ja,                                                                 
!   ia  = input matrix in CSR format                                    
! perm         = integer array of length n containing the indices of the
!         to be extracted.                                              
!                                                                       
! job   = job indicator. if (job  /= 1) values are not copied (i.e.,    
!         only pattern is copied).                                      
!                                                                       
! on return:                                                            
!-----------                                                            
! b,ja,                                                                 
! ib   = matrix in csr format. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1)   
!     contain the value and column indices respectively of the nnz      
!     nonzero elements of the permuted matrix. thus ib(1)=ipos.         
!                                                                       
! Notes:                                                                
!-------                                                                
!  algorithm is NOT in place                                            
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
      integer ko, irow, k, i 
      logical values 
!-----------------------------------------------------------------------
      values = (job  ==  1) 
      ko = ipos 
      ib(1) = ko 
      do 900 i=i1,i2 
         irow = perm(i) 
         do 800 k=ia(irow),ia(irow+1)-1 
            if (values) b(ko) = a(k) 
            jb(ko) = ja(k) 
            ko=ko+1 
  800    continue 
         ib(i-i1+2) = ko 
  900 continue 
      return 
      end subroutine  dperm1

!-----------------------------------------------------------------------

      subroutine dperm2 (i1,i2,a,ja,ia,b,jb,ib,cperm,rperm,istart,      &
     &        ipos,job)                                                 
     integer ipos 
     integer i1,i2,job,istart,ja(*),ia(*),jb(*),ib(*),cperm(*),rperm(*) 
      real(wp) a(*),b(*) 
!-----------------------------------------------------------------------
!     general submatrix permutation/ extraction routine.                
!-----------------------------------------------------------------------
!     extracts rows rperm(i1), rperm(i1+1), ..., rperm(i2) and does an  
!     associated column permutation (using array cperm). The resulting  
!     submatrix is constructed in b, jb, ib. For added flexibility, the 
!     extracted elements are put in sequence starting from row 'istart' 
!     of B. In addition a pointer ipos to the beginning of arrays b,jb, 
!     is also allowed (i.e., nonzero elements are accumulated starting i
!     position ipos of b, jb). In most applications istart and ipos are 
!     equal to one. However, the generality adds substantial flexiblity.
!     EXPLE: (1) to permute msr to msr (excluding diagonals)            
!     call dperm2 (1,n,a,ja,ja,b,jb,jb,rperm,rperm,1,n+2)               
!            (2) To extract rows 1 to 10: define rperm and cperm to be  
!     identity permutations (rperm(i)=i, i=1,n) and then                
!            call dperm2 (1,10,a,ja,ia,b,jb,ib,rperm,rperm,1,1)         
!            (3) to achieve a symmetric permutation as defined by perm: 
!            call dperm2 (1,10,a,ja,ia,b,jb,ib,perm,perm,1,1)           
!            (4) to get a symmetric permutation of A and append the     
!            resulting data structure to A's data structure (useful!)   
!            call dperm2 (1,10,a,ja,ia,a,ja,ia(n+1),perm,perm,1,ia(n+1))
!-----------------------------------------------------------------------
! Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991.                          
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n         = dimension of the matrix                                   
! i1,i2 = extract rows rperm(i1) to rperm(i2) of A, with i1<i2.         
!                                                                       
! a,ja,                                                                 
!   ia  = input matrix in CSR format                                    
! cperm = integer array of length n containing the permutation arrays   
!          for the columns: cperm(i) is the destination of column j,    
!         i.e., any column index ja(k) is transformed into cperm(ja(k)) 
!                                                                       
! rperm        =  permutation array for the rows. rperm(i) = origin (in 
!          row i in B. This is the reverse permutation relative to the  
!          ones used in routines cperm, dperm,....                      
!          rows rperm(i1), rperm(i1)+1, ... rperm(i2) are               
!          extracted from A and stacked into B, starting in row istart  
!          of B.                                                        
! istart= starting row for B where extracted matrix is to be added.     
!         this is also only a pointer of the be beginning address for   
!         ib , on return.                                               
! ipos  = beginning position in arrays b and jb where to start copying  
!         elements. Thus, ib(istart) = ipos.                            
!                                                                       
! job   = job indicator. if (job  /= 1) values are not copied (i.e.,    
!         only pattern is copied).                                      
!                                                                       
! on return:                                                            
!-----------                                                            
! b,ja,                                                                 
! ib   = matrix in csr format. positions 1,2,...,istart-1 of ib         
!     are not touched. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1)           
!     contain the value and column indices respectively of the nnz      
!     nonzero elements of the permuted matrix. thus ib(istart)=ipos.    
!                                                                       
! Notes:                                                                
!-------                                                                
!  1) algorithm is NOT in place                                         
!  2) column indices may not be sorted on return even  though they      
!     may be on entry.                                                  
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
      integer ko, irow, k, i 
      logical values 
!-----------------------------------------------------------------------
      values = (job  ==  1) 
      ko = ipos 
      ib(istart) = ko 
      do 900 i=i1,i2 
         irow = rperm(i) 
         do 800 k=ia(irow),ia(irow+1)-1 
            if (values) b(ko) = a(k) 
            jb(ko) = cperm(ja(k)) 
            ko=ko+1 
  800    continue 
         ib(istart+i-i1+1) = ko 
  900 continue 
      return 
      end subroutine  dperm2

!-----------------------------------------------------------------------

      subroutine dmperm (nrow,a,ja,ao,jao,perm,job) 
      integer nrow,ja(*),jao(*),perm(nrow),job 
      real(wp) a(*),ao(*) 
!-----------------------------------------------------------------------
! This routine performs a symmetric permutation of the rows and         
! columns of a matrix stored in MSR format. i.e., it computes           
! B = P A transp(P), where P, is  a permutation matrix.                 
! P maps row i into row perm(i) and column j into column perm(j):       
!      a(i,j)    becomes   a(perm(i),perm(j)) in new matrix             
! (i.e.  ao(perm(i),perm(j)) = a(i,j) )                                 
! calls dperm.                                                          
!-----------------------------------------------------------------------
! Y. Saad, Nov 15, 1991.                                                
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n         = dimension of the matrix                                   
! a, ja = input matrix in MSR format.                                   
! perm         = integer array of length n containing the permutation ar
!          for the rows: perm(i) is the destination of row i in the     
!         permuted matrix -- also the destination of column i in case   
!         permutation is symmetric (job  <=  2)                         
!                                                                       
! job        = integer indicating the work to be done:                  
!                 job = 1        permute a, ja, ia into ao, jao, iao    
!                 job = 2 permute matrix ignoring real values.          
!                                                                       
! on return:                                                            
!-----------                                                            
! ao, jao = output matrix in MSR.                                       
!                                                                       
! in case job  ==  2 a and ao are never referred to and can be dummy    
! arguments.                                                            
!                                                                       
! Notes:                                                                
!-------                                                                
!  1) algorithm is NOT in place                                         
!  2) column indices may not be sorted on return even  though they may b
!     on entry.                                                         
!----------------------------------------------------------------------c
!     local variables                                                   
!                                                                       
      integer n1, n2, j 
      n1 = nrow+1 
      n2 = n1+1 
!                                                                       
      call dperm (nrow,a,ja,ja,ao(n2),jao(n2),jao,perm,perm,job) 
!                                                                       
      jao(1) = n2 
      do 101 j=1, nrow 
         ao(perm(j)) = a(j) 
         jao(j+1) = jao(j+1)+n1 
  101 continue 
      return
      end subroutine dmperm

!-----------------------------------------------------------------------
  subroutine dvperm ( n, x, perm )

  !*****************************************************************************80
  !
  !! DVPERM performs an in-place permutation of a real vector.
  !
  !  Discussion:
  !
  !    This routine permutes a real vector X using a permutation PERM.
  !
  !    On return, the vector X satisfies,
  !
  !      x(perm(j)) :== x(j), j = 1,2,.., n
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the length of X.
  !
  !    Input/output, real X(N), the vector to be permuted.
  !
  !    Input, integer ( kind = 4 ) PERM(N), the permutation.
  !
    implicit none

    integer :: n

    integer :: ii
    integer :: init
    integer :: k
    integer :: next
    integer :: perm(n)
    real(wp) :: tmp
    real(wp) :: tmp1
    real(wp) :: x(n)

    init = 1
    tmp = x(init)
    ii = perm(init)
    perm(init)= -perm(init)
    k = 0
  !
  !  The main loop.
  !
   6  continue

     k = k + 1
  !
  !  Save the chased element.
  !
    tmp1 = x(ii)
    x(ii) = tmp
    next = perm(ii)

    if ( next < 0 ) then
      go to 65
    end if
  !
  !  Test for end.
  !
    if ( n < k ) then
      perm(1:n) = -perm(1:n)
      return
    end if

    tmp = tmp1
    perm(ii) = -perm(ii)
    ii = next
  !
  !  End of the loop.
  !
    go to 6
  !
  !  Reinitialize cycle.
  !
   65   continue

    init = init + 1

    if ( n < init ) then 
      perm(1:n) = -perm(1:n)
      return
    end if

    if ( perm(init) < 0 ) then
      go to 65
    end if

    tmp = x(init)
    ii = perm(init)
    perm(init) = -perm(init)
    go to 6

  end subroutine dvperm

!-----------------------------------------------------------------------

      subroutine ivperm (n, ix, perm) 
      integer n, perm(n), ix(n) 
!-----------------------------------------------------------------------
! this subroutine performs an in-place permutation of an integer vector 
! ix according to the permutation array perm(*), i.e., on return,       
! the vector x satisfies,                                               
!                                                                       
!        ix(perm(j)) :== ix(j), j=1,2,.., n                             
!                                                                       
!-----------------------------------------------------------------------
! on entry:                                                             
!---------                                                              
! n         = length of vector x.                                       
! perm         = integer array of length n containing the permutation  a
! ix        = input vector                                              
!                                                                       
! on return:                                                            
!----------                                                             
! ix        = vector x permuted according to ix(perm(*)) :=  ix(*)      
!                                                                       
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables                                                       
      integer tmp, tmp1, ii, init, j, k, next 
!                                                                       
      init      = 1 
      tmp        = ix(init) 
      ii        = perm(init) 
      perm(init)= -perm(init) 
      k         = 0 
!                                                                       
! loop                                                                  
!                                                                       
    6 k = k+1 
!                                                                       
! save the chased element --                                            
!                                                                       
      tmp1          = ix(ii) 
      ix(ii)     = tmp 
      next          = perm(ii) 
      if (next  <  0 ) goto 65 
!                                                                       
! test for end                                                          
!                                                                       
      if (k  >  n) goto 101 
      tmp       = tmp1 
      perm(ii)  = - perm(ii) 
      ii        = next 
!                                                                       
! end loop                                                              
!                                                                       
      goto 6 
!                                                                       
! reinitilaize cycle --                                                 
!                                                                       
   65 init      = init+1 
      if (init  >  n) goto 101 
      if (perm(init)  <  0) goto 65 
      tmp        = ix(init) 
      ii        = perm(init) 
      perm(init)=-perm(init) 
      goto 6 
!                                                                       
  101 continue 
      do 200 j=1, n 
         perm(j) = -perm(j) 
  200 continue 
!                                                                       
      return 
      end subroutine  ivperm

!-----------------------------------------------------------------------

      subroutine retmx (n,a,ja,ia,dd) 
      real(wp) a(*),dd(*) 
      integer n,ia(*),ja(*) 
!-----------------------------------------------------------------------
! returns in dd(*) the max absolute value of elements in row *.         
! used for scaling purposes. superseded by rnrms  .                     
!                                                                       
! on entry:                                                             
! n        = dimension of A                                             
! a,ja,ia                                                               
!        = matrix stored in compressed sparse row format                
! dd        = real(wp) array of length n. On output,entry dd(i) contains  
!          the element of row i that has the largest absolute value.    
!          Moreover the sign of dd is modified such that it is the      
!          same as that of the diagonal element in row i.               
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
! local variables                                                       
      integer k2, i, k1, k 
      real(wp) t, t1, t2

      t = 0.0_wp
      t1 = 0.0_wp
      t2 = 0.0_wp
!                                                                       
! initialize                                                            
!                                                                       
      k2 = 1 
      do 11 i=1,n 
         k1 = k2 
         k2 = ia(i+1) - 1 
         t = 0.0d0 
         do 101  k=k1,k2 
            t1 = abs(a(k)) 
            if (t1  >  t) t = t1 
            if (ja(k)  ==  i) then 
               if (a(k)  >=  0.0) then 
                  t2 = a(k) 
               else 
                  t2 = - a(k) 
               endif 
            endif 
  101    continue 
         dd(i) =  t2*t 
!     we do not invert diag                                             
   11 continue 
      return 
      end subroutine  retmx

!-----------------------------------------------------------------------

      subroutine diapos  (n,ja,ia,idiag) 
      integer n
      integer ia(n+1), ja(*), idiag(n) 
!-----------------------------------------------------------------------
! this subroutine returns the positions of the diagonal elements of a   
! sparse matrix a, ja, ia, in the array idiag.                          
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
!                                                                       
! n        = integer. row dimension of the matrix a.                    
! a,ja,                                                                 
!    ia = matrix stored compressed sparse row format. a array skipped.  
!                                                                       
! on return:                                                            
!-----------                                                            
! idiag  = integer array of length n. The i-th entry of idiag           
!          points to the diagonal element a(i,i) in the arrays          
!          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)      
!          if no diagonal element is found the entry is set to 0.       
!----------------------------------------------------------------------c
!           Y. Saad, March, 1990                                        
!----------------------------------------------------------------------c
      integer i, k

      do 1 i=1, n 
         idiag(i) = 0 
    1 continue 
!                                                                       
!     sweep through data structure.                                     
!                                                                       
      do  6 i=1,n 
         do 51 k= ia(i),ia(i+1) -1 
            if (ja(k)  ==  i) idiag(i) = k 
   51    continue 
    6 continue 
      return 
      end subroutine  diapos

!-----------------------------------------------------------------------

      subroutine dscaldg (n,a,ja,ia,diag,job) 
      integer n, job
      real(wp) a(*), diag(*),t 
      integer ia(*),ja(*) 
!-----------------------------------------------------------------------
! scales rows by diag where diag is either given (job=0)                
! or to be computed:                                                    
!  job = 1 ,scale row i by  by  +/- max |a(i,j) | and put inverse of    
!       scaling factor in diag(i),where +/- is the sign of a(i,i).      
!  job = 2 scale by 2-norm of each row..                                
! if diag(i) = 0,then diag(i) is replaced by one                        
! (no scaling)..                                                        
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      integer i, j, k, k1, k2

      goto (12,11,10) job+1 
   10 do 110 j=1,n 
         k1= ia(j) 
         k2 = ia(j+1)-1 
         t = 0.0d0 
         do 111 k = k1,k2 
  111       t = t+a(k)*a(k) 
  110       diag(j) = sqrt(t) 
            goto 12 
   11 continue 
      call retmx (n,a,ja,ia,diag) 
!------                                                                 
   12 do 1 j=1,n 
         if (diag(j)  /=  0.0d0) then 
            diag(j) = 1.0d0/diag(j) 
         else 
            diag(j) = 1.0d0 
         endif 
    1 continue 
      do 2 i=1,n 
         t = diag(i) 
         do 21 k=ia(i),ia(i+1) -1 
            a(k) = a(k)*t 
   21    continue 
    2 continue 
      return 
      end subroutine  dscaldg

!-----------------------------------------------------------------------

      subroutine extbdg (n,a,ja,ia,bdiag,nblk,ao,jao,iao) 
      implicit real(wp) (a-h,o-z) 
      integer n, nblk
      real(wp) bdiag(*),a(*),ao(*) 
      integer ia(*),ja(*),jao(*),iao(*) 
!-----------------------------------------------------------------------
! this subroutine extracts the main diagonal blocks of a                
! matrix stored in compressed sparse row format and puts the result     
! into the array bdiag and the remainder in ao,jao,iao.                 
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n        = integer. The row dimension of the matrix a.                
! a,                                                                    
! ja,                                                                   
! ia    = matrix stored in csr format                                   
! nblk  = dimension of each diagonal block. The diagonal blocks are     
!         stored in compressed format rowwise,i.e.,we store in          
!          succession the i nonzeros of the i-th row after those of     
!          row number i-1..                                             
!                                                                       
! on return:                                                            
!----------                                                             
! bdiag = real(wp) array of size (n x nblk) containing the diagonal       
!          blocks of A on return                                        
! ao,                                                                   
! jao,                                                                  
! iao   = remainder of the matrix stored in csr format.                 
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      integer i, j, j1, j2, jj, k, kb, ko, l, ltr, m 

      m = 1 + (n-1)/nblk 
! this version is sequential -- there is a more parallel version        
! that goes through the structure twice ....                            
      ltr =  ((nblk-1)*nblk)/2 
      l = m * ltr 
      do 1 i=1,l 
         bdiag(i) = 0.0d0 
    1 continue 
      ko = 0 
      kb = 1 
      iao(1) = 1 
!-------------------------                                              
      do 11 jj = 1,m 
         j1 = (jj-1)*nblk+1 
         j2 =  min0 (n,j1+nblk-1) 
         do 12 j=j1,j2 
            do 13 i=ia(j),ia(j+1) -1 
               k = ja(i) 
               if (k  <  j1) then 
                  ko = ko+1 
                  ao(ko) = a(i) 
                  jao(ko) = k 
               else if (k  <  j) then 
!     kb = (jj-1)*ltr+((j-j1)*(j-j1-1))/2+k-j1+1                        
!     bdiag(kb) = a(i)                                                  
                  bdiag(kb+k-j1) = a(i) 
               endif 
   13       continue 
            kb = kb + j-j1 
            iao(j+1) = ko+1 
   12    continue 
   11 continue 
      return 
      end subroutine  extbdg

!-----------------------------------------------------------------------

      subroutine getbwd(n,a,ja,ia,ml,mu) 
!-----------------------------------------------------------------------
! gets the bandwidth of lower part and upper part of A.                 
! does not assume that A is sorted.                                     
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n        = integer = the row dimension of the matrix                  
! a, ja,                                                                
!    ia = matrix in compressed sparse row format.                       
!                                                                       
! on return:                                                            
!-----------                                                            
! ml        = integer. The bandwidth of the strict lower part of A      
! mu        = integer. The bandwidth of the strict upper part of A      
!                                                                       
! Notes:                                                                
! ===== ml and mu are allowed to be negative or return. This may be     
!       useful since it will tell us whether a band is confined         
!       in the strict  upper/lower triangular part.                     
!       indeed the definitions of ml and mu are                         
!                                                                       
!       ml = max ( (i-j)  s.t. a(i,j)  /=  0  )                         
!       mu = max ( (j-i)  s.t. a(i,j)  /=  0  )                         
!----------------------------------------------------------------------c
! Y. Saad, Sep. 21 1989                                                c
!----------------------------------------------------------------------c
      integer n 
      real(wp) a(*) 
      integer ja(*),ia(n+1),ml,mu,ldist,i,k 
      ml = - n 
      mu = - n 
      do 3 i=1,n 
         do 31 k=ia(i),ia(i+1)-1 
            ldist = i-ja(k) 
            ml = max(ml,ldist) 
            mu = max(mu,-ldist) 
   31    continue 
    3 continue 
      return 
      end subroutine  getbwd

!-----------------------------------------------------------------------

      subroutine blkfnd (nrow,ja,ia,nblk) 
!-----------------------------------------------------------------------
! This routine attemptps to determine whether or not  the input         
! matrix has a block structure and finds the blocks size                
! if it does. A block matrix is one which is                            
! comprised of small square dense blocks. If there are zero             
! elements within the square blocks and the original data structure     
! takes these zeros into account then blkchk may fail to find the       
! correct block size.                                                   
!-----------------------------------------------------------------------
! on entry                                                              
!---------                                                              
! nrow        = integer equal to the row dimension of the matrix.       
! ja    = integer array containing the column indices of the entries    
!         nonzero entries of the matrix stored by row.                  
! ia    = integer array of length nrow + 1 containing the pointers      
!         beginning of each row in array ja.                            
!                                                                       
! nblk  = integer containing the assumed value of nblk if job = 0       
!                                                                       
! on return                                                             
!----------                                                             
! nblk  = integer containing the value found for nblk when job = 1.     
!         if imsg  /=  0 this value is meaningless however.             
!                                                                       
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      integer nrow, nblk
      integer ia(nrow+1),ja(*) 
!-----------------------------------------------------------------------
! first part of code will find candidate block sizes.                   
! criterion used here is a simple one: scan rows and  determine groups  
! of rows that have the same length and such that the first column      
! number and the last column number are identical.                      
!-----------------------------------------------------------------------
      integer i, i1, i2, iblk, imsg, irow, jf, jl, jlast, jfirst, jrow, len, &
            & len0, minlen

      minlen = ia(2)-ia(1) 
      irow   = 1 
      do 1 i=2,nrow 
         len = ia(i+1)-ia(i) 
         if (len  <  minlen) then 
            minlen = len 
            irow = i 
         endif 
    1 continue 
!                                                                       
!     ---- candidates are all dividers of minlen                        
!                                                                       
      nblk = 1 
      if (minlen  <=  1) return 
!                                                                       
      do 99 iblk = minlen, 1, -1 
         if (mod(minlen,iblk)  /=  0) goto 99 
         len = ia(2) - ia(1) 
         len0 = len 
         jfirst = ja(1) 
         jlast = ja(ia(2)-1) 
         do 10 jrow = irow+1,irow+nblk-1 
            i1 = ia(jrow) 
            i2 = ia(jrow+1)-1 
            len = i2+1-i1 
            jf = ja(i1) 
            jl = ja(i2) 
            if (len  /=  len0 .or. jf  /=  jfirst .or.                  &
     &           jl  /=  jlast) goto 99                                 
   10    continue 
!                                                                       
!     check for this candidate ----                                     
!                                                                       
         call blkchk (nrow,ja,ia,iblk,imsg) 
         if (imsg  ==  0) then 
!                                                                       
!     block size found                                                  
!                                                                       
            nblk = iblk 
            return 
         endif 
   99 continue 
      end subroutine  blkfnd

!-----------------------------------------------------------------------

      subroutine blkchk (nrow,ja,ia,nblk,imsg) 
!-----------------------------------------------------------------------
! This routine checks whether the input matrix is a block               
! matrix with block size of nblk. A block matrix is one which is        
! comprised of small square dense blocks. If there are zero             
! elements within the square blocks and the data structure              
! takes them into account then blkchk may fail to find the              
! correct block size.                                                   
!-----------------------------------------------------------------------
! on entry                                                              
!---------                                                              
! nrow        = integer equal to the row dimension of the matrix.       
! ja    = integer array containing the column indices of the entries    
!         nonzero entries of the matrix stored by row.                  
! ia    = integer array of length nrow + 1 containing the pointers      
!         beginning of each row in array ja.                            
!                                                                       
! nblk  = integer containing the value of nblk to be checked.           
!                                                                       
! on return                                                             
!----------                                                             
!                                                                       
! imsg  = integer containing a message  with the following meaning.     
!          imsg = 0 means that the output value of nblk is a correct    
!                   block size. nblk  <  0 means nblk not correct      
!                   block size.                                         
!          imsg = -1 : nblk does not divide nrow                        
!          imsg = -2 : a starting element in a row is at wrong position 
!             (j  /=  mult*nblk +1 )                                    
!          imsg = -3 : nblk does divide a row length -                  
!          imsg = -4 : an element is isolated outside a block or        
!             two rows in same group have different lengths             
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      integer nrow, nblk, imsg
      integer ia(nrow+1),ja(*) 
!---------------------------------------------------------------------- 
! first part of code will find candidate block sizes.                   
! this is not guaranteed to work . so a check is done at the end        
! the criterion used here is a simple one:                              
! scan rows and determine groups of rows that have the same length      
! and such that the first column number and the last column number      
! are identical.                                                        
!---------------------------------------------------------------------- 
      integer i, i1, ii, irow, j, j2, jstart, k, len, lena, nr

      imsg = 0 
      if (nblk  <=  1) return 
      nr = nrow/nblk 
      if (nr*nblk  /=  nrow) goto 101 
!--   main loop --------------------------------------------------------
      irow = 1 
      do 20 ii=1, nr 
!     i1= starting position for group of nblk rows in original matrix   
         i1 = ia(irow) 
         j2 = i1 
!     lena = length of each row in that group  in the original matrix   
         lena = ia(irow+1)-i1 
!     len = length of each block-row in that group in the output matrix 
         len = lena/nblk 
         if (len* nblk  /=  lena) goto 103 
!                                                                       
!     for each row                                                      
!                                                                       
         do 6 i = 1, nblk 
            irow = irow + 1 
            if (ia(irow)-ia(irow-1)  /=  lena ) goto 104 
!                                                                       
!     for each block                                                    
!                                                                       
            do 7 k=0, len-1 
               jstart = ja(i1+nblk*k)-1 
               if ( (jstart/nblk)*nblk  /=  jstart) goto 102 
!                                                                       
!     for each column                                                   
!                                                                       
               do 5 j=1, nblk 
                  if (jstart+j  /=  ja(j2) )  goto 104 
                  j2 = j2+1 
    5          continue 
    7       continue 
    6    continue 
   20 continue 
!     went through all loops successfully:                              
      return 
  101 imsg = -1 
      return 
  102 imsg = -2 
      return 
  103 imsg = -3 
      return 
  104 imsg = -4 
      return 
      end subroutine  blkchk

!-----------------------------------------------------------------------

      subroutine infdia (n,ja,ia,ind,idiag) 
      integer n, idiag
      integer ia(*), ind(*), ja(*) 
!-----------------------------------------------------------------------
!     obtains information on the diagonals of A.                        
!-----------------------------------------------------------------------
! this subroutine finds the lengths of each of the 2*n-1 diagonals of A 
! it also outputs the number of nonzero diagonals found.                
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! n        = dimension of the matrix a.                                 
!                                                                       
! a,    ..... not needed here.                                          
! ja,                                                                   
! ia    = matrix stored in csr format                                   
!                                                                       
! on return:                                                            
!-----------                                                            
!                                                                       
! idiag = integer. number of nonzero diagonals found.                   
!                                                                       
! ind   = integer array of length at least 2*n-1. The k-th entry in     
!         ind contains the number of nonzero elements in the diagonal   
!         number k, the numbering beeing from the lowermost diagonal    
!         (bottom-left). In other words ind(k) = length of diagonal     
!         whose offset wrt the main diagonal is = - n + k.              
!----------------------------------------------------------------------c
!           Y. Saad, Sep. 21 1989                                      c
!----------------------------------------------------------------------c
      integer i, j, k, n2      

      n2= n+n-1 
      do 1 i=1,n2 
         ind(i) = 0 
    1 continue 
      do 3 i=1, n 
         do 2 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            ind(n+j-i) = ind(n+j-i) +1 
    2    continue 
    3 continue 
!     count the nonzero ones.                                           
      idiag = 0 
      do 41 k=1, n2 
         if (ind(k)  /=  0) idiag = idiag+1 
   41 continue 
      return 
      end subroutine  infdia

!-----------------------------------------------------------------------

  subroutine amubdg ( nrow, ncol, ncolb, ja, ia, jb, ib, ndegr, nnz, iw )

  !*****************************************************************************80
  !
  !! AMUBDG gets the number of nonzero elements in each row of A * B.
  !
  !  Discussion:
  !
  !    The routine also computes the total number of nonzero elements in A * B.
  !
  !    Method: A' * A = sum [over i = 1, nrow]  a(i)^T a(i)
  !    where a(i) = i-th row of  A.  We must be careful not to add  the
  !    elements already accounted for.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix A.
  !
  !    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix A,
  !    (and the row dimension of B).
  !
  !    Input, integer ( kind = 4 ) NCOLB, the column dimension of the matrix B.
  !
  !    Input, ja, ia= row structure of input matrix A: ja = column indices of
  !    the nonzero elements of A stored by rows.
  !    ia = pointer to beginning of each row in ja.
  !
  !    Input, jb, ib, the row structure of input matrix B: jb = column indices of
  !    the nonzero elements of A stored by rows.
  !    ib is a pointer to beginning of each row in jb.
  !
  !    Output, integer ( kind = 4 ) NDEGR(NROW), contains the degrees (the number
  !    of nonzeros in each row of the matrix A * B.
  !
  !    Output, integer ( kind = 4 ) NNZ, the number of nonzero elements 
  !    found in A * B.
  !
  !    Workspace, integer ( kind = 4 ) IW(NCOLB).
  !
    implicit none

    integer :: ncol
    integer :: ncolb
    integer :: nrow

    integer :: ia(nrow+1)
    integer :: ib(ncol+1)
    integer :: ii
    integer :: iw(ncolb)
    integer :: j
    integer :: ja(*)
    integer :: jb(*)
    integer :: jc
    integer :: jr
    integer :: k
    integer :: last
    integer :: ldg
    integer :: ndegr(nrow)
    integer :: nnz

    iw(1:ncolb) = 0
    ndegr(1:nrow) = 0

    do ii = 1, nrow
  !
  !  For each row of A.
  !
      ldg = 0
  !
  !  End-of-linked list.
  !
      last = -1

      do j = ia(ii), ia(ii+1)-1
  !
  !  Row number to be added.
  !
          jr = ja(j)

          do k = ib(jr), ib(jr+1)-1
             jc = jb(k)
  !
  !  Add one element to the linked list.
  !
             if ( iw(jc) == 0 ) then
                ldg = ldg + 1
                iw(jc) = last
                last = jc
             end if

           end do

      end do

      ndegr(ii) = ldg
  !
  !  Reset IW to zero.
  !
      do k = 1, ldg
        j = iw(last)
        iw(last) = 0
        last = j
       end do

    end do

    nnz = sum ( ndegr(1:nrow) )

    return
  end subroutine amubdg

!-----------------------------------------------------------------------
  subroutine aplbdg ( nrow, ncol, ja, ia, jb, ib, ndegr, nnz, iw )

  !*****************************************************************************80
  !
  !! APLBDG gets the number of nonzero elements in each row of A + B.
  !
  !  Discussion:
  !
  !    It also reports the total number of nonzero elements in A + B.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
  !
  !    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  !    Input, b, jb, ib, matrix B in compressed sparse row format.
  !
  !    Output, integer ( kind = 4 ) NDEGR(NROW), the number of nonzeros in each row 
  !    of the matrix A + B.
  !
  !    Output, integer ( kind = 4 ) NNZ, the total number of nonzero elements found 
  !    in A + B.
  !
  !    Workspace, integer ( kind = 4 ) IW(NCOL).
  !
    implicit none

    integer :: ncol
    integer :: nrow

    integer :: ia(nrow+1)
    integer :: ib(nrow+1)
    integer :: ii
    integer :: iw(ncol)
    integer :: j
    integer :: ja(*)
    integer :: jb(*)
    integer :: jc
    integer :: jr
    integer :: k
    integer :: last
    integer :: ldg
    integer :: ndegr(nrow)
    integer :: nnz

    iw(1:ncol) = 0

    ndegr(1:nrow) = 0

    do ii = 1, nrow

       ldg = 0
  !
  !  End-of-linked list.
  !
       last = -1
  !
  !  Row of A.
  !
       do j = ia(ii), ia(ii+1)-1
          jr = ja(j)
  !
  !  Add element to the linked list.
  !
          ldg = ldg + 1
          iw(jr) = last
          last = jr
       end do
  !
  !  Row of B.
  !
       do j = ib(ii), ib(ii+1)-1

          jc = jb(j)
  !
  !  Add one element to the linked list.
  !
          if ( iw(jc) == 0 ) then
             ldg = ldg + 1
             iw(jc) = last
             last = jc
          end if

       end do
  !
  !  Done with row II.
  !
       ndegr(ii) = ldg
  !
  !  Reset IW to zero.
  !
       do k = 1, ldg
          j = iw(last)
          iw(last) = 0
          last = j
       end do

    end do

    nnz = sum ( ndegr(1:nrow) )

    return
  end subroutine aplbdg

!-----------------------------------------------------------------------

      subroutine rnrms   (nrow, nrm, a, ja, ia, diag) 
      integer nrow, nrm
      real(wp) a(*), diag(nrow), scal 
      integer ja(*), ia(nrow+1) 
!-----------------------------------------------------------------------
! gets the norms of each row of A. (choice of three norms)              
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2        
!                  means the 2-nrm, nrm = 0 means max norm              
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! diag = real vector of length nrow containing the norms                
!                                                                       
!-----------------------------------------------------------------      
      integer ii, k, k1, k2      

      do 1 ii=1,nrow 
!                                                                       
!     compute the norm if each element.                                 
!                                                                       
         scal = 0.0d0 
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         if (nrm  ==  0) then 
            do 2 k=k1, k2 
               scal = max(scal,abs(a(k) ) ) 
    2       continue 
         elseif (nrm  ==  1) then 
            do 3 k=k1, k2 
               scal = scal + abs(a(k) ) 
    3       continue 
         else 
            do 4 k=k1, k2 
               scal = scal+a(k)**2 
    4       continue 
         endif 
         if (nrm  ==  2) scal = sqrt(scal) 
         diag(ii) = scal 
    1 continue 
      return 
      end subroutine  rnrms

!-----------------------------------------------------------------------

      subroutine cnrms   (nrow, nrm, a, ja, ia, diag) 
      integer nrow, nrm
      real(wp) a(*), diag(nrow) 
      integer ja(*), ia(nrow+1) 
!-----------------------------------------------------------------------
! gets the norms of each column of A. (choice of three norms)           
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2        
!                  means the 2-nrm, nrm = 0 means max norm              
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! diag = real vector of length nrow containing the norms                
!                                                                       
!-----------------------------------------------------------------      
      integer ii, j, k, k1, k2  

      do 10 k=1, nrow 
         diag(k) = 0.0d0 
   10 continue 
      do 1 ii=1,nrow 
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         do 2 k=k1, k2 
            j = ja(k) 
!     update the norm of each column                                    
            if (nrm  ==  0) then 
               diag(j) = max(diag(j),abs(a(k) ) ) 
            elseif (nrm  ==  1) then 
               diag(j) = diag(j) + abs(a(k) ) 
            else 
               diag(j) = diag(j)+a(k)**2 
            endif 
    2    continue 
    1 continue 
      if (nrm  /=  2) return 
      do 3 k=1, nrow 
         diag(k) = sqrt(diag(k)) 
    3 continue 
      return 
      end subroutine  cnrms

!-----------------------------------------------------------------------

      subroutine roscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      integer nrow 
      real(wp) a(*), b(*), diag(nrow) 
      integer job,nrm,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
!-----------------------------------------------------------------------
! scales the rows of A such that their norms are one on return          
! 3 choices of norms: 1-norm, 2-norm, max-norm.                         
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! job   = integer. job indicator. Job=0 means get array b only          
!         job = 1 means get b, and the integer arrays ib, jb.           
!                                                                       
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2        
!                  means the 2-nrm, nrm = 0 means max norm              
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! diag = diagonal matrix stored as a vector containing the matrix       
!        by which the rows have been scaled, i.e., on return            
!        we have B = Diag*A.                                            
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        = resulting matrix B in compressed sparse row sparse format.
!                                                                       
! ierr  = error message. ierr=0     : Normal return                     
!                        ierr=i > 0 : Row number i is a zero row.       
! Notes:                                                                
!-------                                                                
! 1)        The column dimension of A is not needed.                    
! 2)        algorithm in place (B can take the place of A).             
!-----------------------------------------------------------------      
      integer j

      call rnrms (nrow,nrm,a,ja,ia,diag) 
      ierr = 0 
      do 1 j=1, nrow 
         if (diag(j)  ==  0.0d0) then 
            ierr = j 
            return 
         else 
            diag(j) = 1.0d0/diag(j) 
         endif 
    1 continue 
      call diamua(nrow,job,a,ja,ia,diag,b,jb,ib) 
      return 
      end subroutine  roscal

!-----------------------------------------------------------------------

      subroutine coscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
!-----------------------------------------------------------------------
      integer nrow, nrm
      real(wp) a(*),b(*),diag(nrow) 
      integer job,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
!-----------------------------------------------------------------------
! scales the columns of A such that their norms are one on return       
! result matrix written on b, or overwritten on A.                      
! 3 choices of norms: 1-norm, 2-norm, max-norm. in place.               
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! job   = integer. job indicator. Job=0 means get array b only          
!         job = 1 means get b, and the integer arrays ib, jb.           
!                                                                       
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2        
!                  means the 2-nrm, nrm = 0 means max norm              
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! diag = diagonal matrix stored as a vector containing the matrix       
!        by which the columns have been scaled, i.e., on return         
!        we have B = A * Diag                                           
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        = resulting matrix B in compressed sparse row sparse format.
!                                                                       
! ierr  = error message. ierr=0     : Normal return                     
!                        ierr=i > 0 : Column number i is a zero row.    
! Notes:                                                                
!-------                                                                
! 1)     The column dimension of A is not needed.                       
! 2)     algorithm in place (B can take the place of A).                
!-----------------------------------------------------------------      
      integer j

      call cnrms (nrow,nrm,a,ja,ia,diag) 
      ierr = 0 
      do 1 j=1, nrow 
         if (diag(j)  ==  0.0) then 
            ierr = j 
            return 
         else 
            diag(j) = 1.0d0/diag(j) 
         endif 
    1 continue 
      call amudia (nrow,job,a,ja,ia,diag,b,jb,ib) 
      return 
      end subroutine  coscal

!-----------------------------------------------------------------------

      subroutine addblk(nrowa, ncola, a, ja, ia, ipos, jpos, job,       &
     & nrowb, ncolb, b, jb, ib, nrowc, ncolc, c, jc, ic, nzmx, ierr)    
!      implicit none                                                    
      integer nrowa, nrowb, nrowc, ncola, ncolb, ncolc, ipos, jpos 
      integer nzmx, ierr, job 
      integer ja(1:*), ia(1:*), jb(1:*), ib(1:*), jc(1:*), ic(1:*) 
      real(wp) a(1:*), b(1:*), c(1:*) 
!-----------------------------------------------------------------------
!     This subroutine adds a matrix B into a submatrix of A whose       
!     (1,1) element is located in the starting position (ipos, jpos).   
!     The resulting matrix is allowed to be larger than A (and B),      
!     and the resulting dimensions nrowc, ncolc will be redefined       
!     accordingly upon return.                                          
!     The input matrices are assumed to be sorted, i.e. in each row     
!     the column indices appear in ascending order in the CSR format.   
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrowa    = number of rows in A.                                       
! bcola    = number of columns in A.                                    
! a,ja,ia  = Matrix A in compressed sparse row format with entries sorte
! nrowb    = number of rows in B.                                       
! ncolb    = number of columns in B.                                    
! b,jb,ib  = Matrix B in compressed sparse row format with entries sorte
!                                                                       
! nzmax           = integer. The  length of the arrays c and jc. addblk 
!            stop if the number of nonzero elements in the matrix C     
!            exceeds nzmax. See ierr.                                   
!                                                                       
! on return:                                                            
!----------                                                             
! nrowc    = number of rows in C.                                       
! ncolc    = number of columns in C.                                    
! c,jc,ic  = resulting matrix C in compressed sparse row sparse format  
!            with entries sorted ascendly in each row.                  
!                                                                       
! ierr           = integer. serving as error message.                   
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that addblk stopped while computing the     
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! Notes:                                                                
!-------                                                                
!     this will not work if any of the two input matrices is not sorted 
!-----------------------------------------------------------------------
      logical values 
      integer i,j1,j2,ka,kb,kc,kamax,kbmax 
      values = (job  /=  0) 
      ierr = 0 
      nrowc = max(nrowa, nrowb+ipos-1) 
      ncolc = max(ncola, ncolb+jpos-1) 
      kc = 1 
      kamax = 0
      kbmax = 0 
      ic(1) = kc 
!                                                                       
      do 10 i=1, nrowc 
         if (i <= nrowa) then 
            ka = ia(i) 
            kamax = ia(i+1)-1 
         else 
            ka = ia(nrowa+1) 
         end if 
         if ((i >= ipos).and.((i-ipos) <= nrowb)) then 
            kb = ib(i-ipos+1) 
            kbmax = ib(i-ipos+2)-1 
         else 
            kb = ib(nrowb+1) 
         end if 
!                                                                       
!     a do-while type loop -- goes through all the elements in a row.   
!                                                                       
   20    continue 
         if (ka  <=  kamax) then 
            j1 = ja(ka) 
         else 
            j1 = ncolc+1 
         endif 
         if (kb  <=  kbmax) then 
            j2 = jb(kb) + jpos - 1 
         else 
            j2 = ncolc+1 
         endif 
!                                                                       
!     if there are more elements to be added.                           
!                                                                       
         if ((ka  <=  kamax .or. kb  <=  kbmax) .and.                   &
     &        (j1  <=  ncolc .or. j2  <=  ncolc)) then                  
!                                                                       
!     three cases                                                       
!                                                                       
            if (j1  ==  j2) then 
               if (values) c(kc) = a(ka)+b(kb) 
               jc(kc) = j1 
               ka = ka+1 
               kb = kb+1 
               kc = kc+1 
            else if (j1  <  j2) then 
               jc(kc) = j1 
               if (values) c(kc) = a(ka) 
               ka = ka+1 
               kc = kc+1 
            else if (j1  >  j2) then 
               jc(kc) = j2 
               if (values) c(kc) = b(kb) 
               kb = kb+1 
               kc = kc+1 
            endif 
            if (kc  >  nzmx) goto 999 
            goto 20 
         end if 
         ic(i+1) = kc 
   10 continue 
      return 
  999 ierr = i 
      return 
      end subroutine  addblk

!-----------------------------------------------------------------------

      subroutine block_addblk(nb, nrowa, ncola, a, ja, ia, ipos, jpos,  &
     & job, nrowb, ncolb, b, jb, ib, nrowc, ncolc, c, jc, ic, nzmx, ierr)    

      implicit none                                                    

      integer,                intent(in)          :: nb
      integer,                intent(in)          :: nrowa, ncola
      integer,                intent(in)          :: nrowb, ncolb
      integer,                intent(in)          :: ipos,  jpos
      integer,                intent(in)          :: job,   nzmx
      integer,                intent(inout)       :: nrowc, ncolc
      integer,                intent(out)         :: ierr

      integer , dimension(*), intent(in)          :: ia, ja
      integer , dimension(*), intent(in)          :: ib, jb
      integer , dimension(*), intent(inout)       :: ic, jc 

      real(wp), dimension(nb,nb,*), intent(in)    :: a, b
      real(wp), dimension(nb,nb,*), intent(inout) :: c

      logical         :: values 
      integer         ::  i,j1,j2,ka,kb,kc,kamax,kbmax 

!-----------------------------------------------------------------------
!     This subroutine adds a matrix B into a submatrix of A whose       
!     (1,1) element is located in the starting position (ipos, jpos).   
!     The resulting matrix is allowed to be larger than A (and B),      
!     and the resulting dimensions nrowc, ncolc will be redefined       
!     accordingly upon return.                                          
!     The input matrices are assumed to be sorted, i.e. in each row     
!     the column indices appear in ascending order in the CSR format.   
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrowa    = number of rows in A.                                       
! bcola    = number of columns in A.                                    
! a,ja,ia  = Matrix A in compressed sparse row format with entries sorte
! nrowb    = number of rows in B.                                       
! ncolb    = number of columns in B.                                    
! b,jb,ib  = Matrix B in compressed sparse row format with entries sorte
!                                                                       
! nzmax           = integer. The  length of the arrays c and jc. addblk 
!            stop if the number of nonzero elements in the matrix C     
!            exceeds nzmax. See ierr.                                   
!                                                                       
! on return:                                                            
!----------                                                             
! nrowc    = number of rows in C.                                       
! ncolc    = number of columns in C.                                    
! c,jc,ic  = resulting matrix C in compressed sparse row sparse format  
!            with entries sorted ascendly in each row.                  
!                                                                       
! ierr           = integer. serving as error message.                   
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that addblk stopped while computing the     
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! Notes:                                                                
!-------                                                                
!     this will not work if any of the two input matrices is not sorted 
!-----------------------------------------------------------------------

      values = (job  /=  0) 
      ierr = 0 
      nrowc = max(nrowa, nrowb+ipos-1) 
      ncolc = max(ncola, ncolb+jpos-1) 
      kc = 1 
      kamax = 0
      kbmax = 0 
      ic(1) = kc 
!                                                                       
      do 10 i=1, nrowc 
         if (i <= nrowa) then 
            ka = ia(i) 
            kamax = ia(i+1)-1 
         else 
            ka = ia(nrowa+1) 
         end if 
         if ((i >= ipos).and.((i-ipos) <= nrowb)) then 
            kb = ib(i-ipos+1) 
            kbmax = ib(i-ipos+2)-1 
         else 
            kb = ib(nrowb+1) 
         end if 
!                                                                       
!     a do-while type loop -- goes through all the elements in a row.   
!                                                                       
   20    continue 
         if (ka  <=  kamax) then 
            j1 = ja(ka) 
         else 
            j1 = ncolc+1 
         endif 
         if (kb  <=  kbmax) then 
            j2 = jb(kb) + jpos - 1 
         else 
            j2 = ncolc+1 
         endif 
!                                                                       
!     if there are more elements to be added.                           
!                                                                       
         if ((ka  <=  kamax .or. kb  <=  kbmax) .and.                   &
     &        (j1  <=  ncolc .or. j2  <=  ncolc)) then                  
!                                                                       
!     three cases                                                       
!                                                                       
            if (j1  ==  j2) then 
               if (values) c(:,:,kc) = a(:,:,ka)+b(:,:,kb) 
               jc(kc) = j1 
               ka = ka+1 
               kb = kb+1 
               kc = kc+1 
            else if (j1  <  j2) then 
               jc(kc) = j1 
               if (values) c(:,:,kc) = a(:,:,ka) 
               ka = ka+1 
               kc = kc+1 
            else if (j1  >  j2) then 
               jc(kc) = j2 
               if (values) c(:,:,kc) = b(:,:,kb) 
               kb = kb+1 
               kc = kc+1 
            endif 
            if (kc  >  nzmx) goto 999 
            goto 20 
         end if 
         ic(i+1) = kc 
   10 continue 
      return 
  999 ierr = i 
      return 
      end subroutine  block_addblk

!-----------------------------------------------------------------------

      subroutine get1up (n,ja,ia,ju) 
      integer  n, ja(*),ia(*),ju(*) 
!---------------------------------------------------------------------- 
! obtains the first element of each row of the upper triangular part    
! of a matrix. Assumes that the matrix is already sorted.               
!-----------------------------------------------------------------------
! parameters                                                            
! input                                                                 
! -----                                                                 
! ja      = integer array containing the column indices of aij          
! ia      = pointer array. ia(j) contains the position of the           
!           beginning of row j in ja                                    
!                                                                       
! output                                                                
! ------                                                                
! ju      = integer array of length n. ju(i) is the address in ja       
!           of the first element of the uper triangular part of         
!           of A (including rthe diagonal. Thus if row i does have      
!           a nonzero diagonal element then ju(i) will point to it.     
!           This is a more general version of diapos.                   
!-----------------------------------------------------------------------
! local vAriables                                                       
      integer i, k 
!                                                                       
      do 5 i=1, n 
         ju(i) = 0 
         k = ia(i) 
!                                                                       
    1    continue 
         if (ja(k)  >=  i) then 
            ju(i) = k 
            goto 5 
         elseif (k  <  ia(i+1) -1) then 
            k=k+1 
!                                                                       
! go try next element in row                                            
!                                                                       
            goto 1 
         endif 
    5 continue 
      return 
      end subroutine  get1up

!---------------------------------------------------------------------- 

      subroutine xtrows (i1,i2,a,ja,ia,ao,jao,iao,iperm,job) 
      integer i1,i2,ja(*),ia(*),jao(*),iao(*),iperm(*),job 
      real(wp) a(*),ao(*) 
!-----------------------------------------------------------------------
! this subroutine extracts given rows from a matrix in CSR format.      
! Specifically, rows number iperm(i1), iperm(i1+1), ...., iperm(i2)     
! are extracted and put in the output matrix ao, jao, iao, in CSR       
! format.  NOT in place.                                                
! Youcef Saad -- coded Feb 15, 1992.                                    
!-----------------------------------------------------------------------
! on entry:                                                             
!----------                                                             
! i1,i2   = two integers indicating the rows to be extracted.           
!           xtrows will extract rows iperm(i1), iperm(i1+1),..,iperm(i2)
!           from original matrix and stack them in output matrix        
!           ao, jao, iao in csr format                                  
!                                                                       
! a, ja, ia = input matrix in csr format                                
!                                                                       
! iperm        = integer array of length nrow containing the reverse per
!         array for the rows. row number iperm(j) in permuted matrix PA 
!         used to be row number j in unpermuted matrix.                 
!         ---> a(i,j) in the permuted matrix was a(iperm(i),j)          
!         in the inout matrix.                                          
!                                                                       
! job        = integer indicating the work to be done:                  
!                 job  /=  1 : get structure only of output matrix,,    
!               i.e., ignore real values. (in which case arrays a       
!               and ao are not used nor accessed).                      
!                 job = 1        get complete data structure of output m
!               (i.e., including arrays ao and iao).                    
!------------                                                           
! on return:                                                            
!------------                                                           
! ao, jao, iao = input matrix in a, ja, ia format                       
! note :                                                                
!        if (job /= 1)  then the arrays a and ao are not used.          
!----------------------------------------------------------------------c
!           Y. Saad, revised May  2, 1990                              c
!----------------------------------------------------------------------c
      logical values 
      integer ii, j, k, ko 

      values = (job  ==  1) 
!                                                                       
! copying                                                               
!                                                                       
      ko = 1 
      iao(1) = ko 
      do 100 j=i1,i2 
!                                                                       
! ii=iperm(j) is the index of old row to be copied.                     
!                                                                       
         ii = iperm(j) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k) 
            ko = ko+1 
   60    continue 
         iao(j-i1+2) = ko 
  100 continue 
!                                                                       
      return 
      end subroutine  xtrows

!-----------------------------------------------------------------------

      subroutine csrkvstr(n, ia, ja, nr, kvstr) 
!-----------------------------------------------------------------------
      integer n, ia(n+1), ja(*), nr, kvstr(*) 
!-----------------------------------------------------------------------
!     Finds block row partitioning of matrix in CSR format.             
!-----------------------------------------------------------------------
!     On entry:                                                         
!--------------                                                         
!     n       = number of matrix scalar rows                            
!     ia,ja   = input matrix sparsity structure in CSR format           
!                                                                       
!     On return:                                                        
!---------------                                                        
!     nr      = number of block rows                                    
!     kvstr   = first row number for each block row                     
!                                                                       
!     Notes:                                                            
!-----------                                                            
!     Assumes that the matrix is sorted by columns.                     
!     This routine does not need any workspace.                         
!                                                                       
!-----------------------------------------------------------------------
!     local variables                                                   
      integer i, j, jdiff 
!-----------------------------------------------------------------------
      nr = 1 
      kvstr(1) = 1 
!---------------------------------                                      
      do i = 2, n 
         jdiff = ia(i+1)-ia(i) 
         if (jdiff  ==  ia(i)-ia(i-1)) then 
            do j = ia(i), ia(i+1)-1 
               if (ja(j)  /=  ja(j-jdiff)) then 
                  nr = nr + 1 
                  kvstr(nr) = i 
                  goto 299 
               endif 
            enddo 
  299       continue 
         else 
! 300       nr = nr + 1    MHC
            nr = nr + 1 
            kvstr(nr) = i 
         endif 
      enddo 
      kvstr(nr+1) = n+1 
!---------------------------------                                      
      return 
      end subroutine  csrkvstr

      subroutine csrkvstc(n, ia, ja, nc, kvstc, iwk) 

!-----------------------------------------------------------------------
      integer n, ia(n+1), ja(*), nc, kvstc(*), iwk(*) 
!-----------------------------------------------------------------------
!     Finds block column partitioning of matrix in CSR format.          
!-----------------------------------------------------------------------
!     On entry:                                                         
!--------------                                                         
!     n       = number of matrix scalar rows                            
!     ia,ja   = input matrix sparsity structure in CSR format           
!                                                                       
!     On return:                                                        
!---------------                                                        
!     nc      = number of block columns                                 
!     kvstc   = first column number for each block column               
!                                                                       
!     Work space:                                                       
!----------------                                                       
!     iwk(*) of size equal to the number of scalar columns plus one.    
!        Assumed initialized to 0, and left initialized on return.      
!                                                                       
!     Notes:                                                            
!-----------                                                            
!     Assumes that the matrix is sorted by columns.                     
!                                                                       
!-----------------------------------------------------------------------
!     local variables                                                   
      integer i, j, k, ncol 
!                                                                       
!-----------------------------------------------------------------------
!-----use ncol to find maximum scalar column number                     
      ncol = 0 
!-----mark the beginning position of the blocks in iwk                  
      do i = 1, n 
         if (ia(i)  <  ia(i+1)) then 
            j = ja(ia(i)) 
            iwk(j) = 1 
            do k = ia(i)+1, ia(i+1)-1 
               j = ja(k) 
               if (ja(k-1) /= j-1) then 
                  iwk(j) = 1 
                  iwk(ja(k-1)+1) = 1 
               endif 
            enddo 
            iwk(j+1) = 1 
            ncol = max0(ncol, j) 
         endif 
      enddo 
!---------------------------------                                      
      nc = 1 
      kvstc(1) = 1 
      do i = 2, ncol+1 
         if (iwk(i) /= 0) then 
            nc = nc + 1 
            kvstc(nc) = i 
            iwk(i) = 0 
         endif 
      enddo 
      nc = nc - 1 
!---------------------------------                                      
      return 
      end subroutine  csrkvstc

!-----------------------------------------------------------------------

      subroutine kvstmerge(nr, kvstr, nc, kvstc, n, kvst) 
!-----------------------------------------------------------------------
      integer nr, kvstr(nr+1), nc, kvstc(nc+1), n, kvst(*) 
!-----------------------------------------------------------------------
!     Merges block partitionings, for conformal row/col pattern.        
!-----------------------------------------------------------------------
!     On entry:                                                         
!--------------                                                         
!     nr,nc   = matrix block row and block column dimension             
!     kvstr   = first row number for each block row                     
!     kvstc   = first column number for each block column               
!                                                                       
!     On return:                                                        
!---------------                                                        
!     n       = conformal row/col matrix block dimension                
!     kvst    = conformal row/col block partitioning                    
!                                                                       
!     Notes:                                                            
!-----------                                                            
!     If matrix is not square, this routine returns without warning.    
!                                                                       
!-----------------------------------------------------------------------
!-----local variables                                                   
      integer i,j 
!---------------------------------                                      
      if (kvstr(nr+1)  /=  kvstc(nc+1)) return 
      i = 1 
      j = 1 
      n = 1 
  200 if (i  >  nr+1) then 
         kvst(n) = kvstc(j) 
         j = j + 1 
      elseif (j  >  nc+1) then 
         kvst(n) = kvstr(i) 
         i = i + 1 
      elseif (kvstc(j)  ==  kvstr(i)) then 
         kvst(n) = kvstc(j) 
         j = j + 1 
         i = i + 1 
      elseif (kvstc(j)  <  kvstr(i)) then 
         kvst(n) = kvstc(j) 
         j = j + 1 
      else 
         kvst(n) = kvstr(i) 
         i = i + 1 
      endif 
      n = n + 1 
      if (i <= nr+1 .or. j <= nc+1) goto 200 
      n = n - 2 
!---------------------------------                                      
      return 
      end subroutine  kvstmerge

!-----------------------------------------------------------------------
  subroutine amub ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
    iw, ierr )

  !*****************************************************************************80
  !
  !! AMUB performs the matrix product C = A * B.
  !
  !  Discussion:
  !
  !    The column dimension of B is not needed.
  !
  !  Modified:
  !
  !    08 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
  !
  !    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
  !
  !    Input, integer ( kind = 4 ) JOB, job indicator.  When JOB = 0, only the
  !    structure is computed, that is, the arrays JC and IC, but the real values
  !    are ignored.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  !    Input, b, jb, ib, matrix B in compressed sparse row format.
  !
  !    Input, integer ( kind = 4 ) NZMAX, the length of the arrays c and jc.
  !    The routine will stop if the result matrix C  has a number
  !    of elements that exceeds exceeds NZMAX.
  !
  ! on return:
  !
  ! c,
  ! jc,
  ! ic    = resulting matrix C in compressed sparse row sparse format.
  !
  ! ierr      = integer ( kind = 4 ). serving as error message.
  !         ierr = 0 means normal return,
  !         ierr > 0 means that amub stopped while computing the
  !         i-th row  of C with i = ierr, because the number
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !
  !  iw      = integer ( kind = 4 ) work array of length equal to the number of
  !         columns in A.
  !
    implicit none

    integer :: ncol
    integer :: nrow
    integer :: nzmax

    real(wp) :: a(*)
    real(wp) :: b(*)
    real(wp) :: c(nzmax)
    integer :: ia(nrow+1)
    integer :: ib(ncol+1)
    integer :: ic(ncol+1)
    integer :: ierr
    integer :: ii
    integer :: iw(ncol)
    integer :: ja(*)
    integer :: jb(*)
    integer :: jc(nzmax)
    integer :: jcol
    integer :: jj
    integer :: job
    integer :: jpos
    integer :: k
    integer :: ka
    integer :: kb
    integer :: len
    real(wp) :: scal
    logical :: values

    values = ( job /= 0 )
    len = 0
    ic(1) = 1
    ierr = 0
    scal = 0.0_wp
  !
  !  Initialize IW.
  !
    iw(1:ncol) = 0

    do ii = 1, nrow
  !
  !  Row I.
  !
      do ka = ia(ii), ia(ii+1)-1

        if ( values ) then
          scal = a(ka)
        end if

        jj = ja(ka)

        do kb = ib(jj), ib(jj+1)-1

             jcol = jb(kb)
             jpos = iw(jcol)

             if ( jpos == 0 ) then
                len = len + 1
                if ( nzmax < len ) then
                   ierr = ii
                   return
                end if
                jc(len) = jcol
                iw(jcol)= len
                if ( values ) then
                  c(len) = scal * b(kb)
                end if
             else
                if ( values ) then
                  c(jpos) = c(jpos) + scal * b(kb)
                end if
             end if

           end do

      end do

      do k = ic(ii), len
        iw(jc(k)) = 0
      end do

      ic(ii+1) = len + 1

    end do

    return
  end subroutine amub

!-----------------------------------------------------------------------

  pure subroutine a_sc_mu_b_bl(nrow,ncol,nb,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)                          

      implicit none

      integer,                      intent(in)    :: nrow, ncol, nb, nzmax
      real(wp), dimension(*),       intent(in)    ::  a
      real(wp), dimension(nb,nb,*), intent(in)    ::  b
      real(wp), dimension(nb,nb,*), intent(out)   ::  c
      
      integer, dimension(nrow+1),   intent(in)    :: ia
      integer, dimension(*),        intent(in)    :: ja
      integer, dimension(*),        intent(in)    :: ib,jb
      integer, dimension(*),        intent(out)   :: ic,jc
 
      integer, dimension(*),       intent(inout) :: iw

      integer,                     intent(out)   :: ierr

      integer                  :: k,ka,jj,ii,j,jcol,jpos,kb,len

      real(wp)                                  :: scal 

!-----------------------------------------------------------------------
! performs the matrix by matrix product C = (A tensor product I_nb) B  and 
! stores the resulting matrix in thr 3D array C whose shape is (nb,nb,nzmax)
!-----------------------------------------------------------------------
! on entry:  
! ---------                                                             
! nrow  = integer. The row dimension of A          
! ncol  = integer. The column dimension of B   
! nb    = integer. The block size of the matrices a, b, and c 
!                                                                       
! a,  ja, ia ; Matrix A in compressed sparse row format.                      
! b,  jb, ib ; Matrix B in compressed sparse row format.                    
!                                                                       
! nzmax = integer. The  length of the arrays c and jc.                  
!         a_sc_mu_b_bl will stop if the result matrix C has a number           
!         of elements that exceeds exceeds nzmax. See ierr.
!                                                                       
! on return: ! ----------                                                             
! c,  jc, ic ; resulting matrix C in compressed sparse row sparse format.  
!              Each element of the vector c is a block of size nb.
!                                                                       
! ierr  = integer. serving as error message.                            
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that amub stopped while computing the       
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! work arrays:                                                          
!------------                                                           
! iw    = integer work array of length equal to the number of columns in A.                                                 

!   The row dimension of B is not needed. However there is no checking  
!   on the condition that ncol(A) = nrow(B).                            

      continue

      len   = 0 
      ic(1) = 1 
      ierr  = 0 

!     initialize array iw.                                              
      do j=1, ncol 
         iw(j) = 0 
      enddo
!                                                                       
      do ii=1, nrow                   !     row i                                                             

         do ka=ia(ii), ia(ii+1)-1 

                 scal =  a(ka) 
                 jj   = ja(ka) 

            do kb=ib(jj),ib(jj+1)-1 
               jcol = jb(kb) 
               jpos = iw(jcol) 
               if (jpos  ==  0) then 
                  len = len+1 
                  if (len  >  nzmax) then 
                     ierr = ii 
                     return 
                  endif 
                  jc(len) = jcol 
                  iw(jcol)= len 
                  c(:,:,len)  = scal*b(:,:,kb)
               else 
                  c(:,:,jpos) = c(:,:,jpos) + scal*b(:,:,kb) 
               endif 
            enddo

         enddo

         do k=ic(ii), len 
            iw(jc(k)) = 0 
         enddo

         ic(ii+1) = len+1 

      enddo

      end subroutine a_sc_mu_b_bl

!-----------------------------------------------------------------------

  pure subroutine a_bl_mu_b_bl(nrow,ncol,nb,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)                          

      implicit none

      integer,                      intent(in)    :: nrow, ncol, nb, nzmax
      real(wp), dimension(nb,nb,*), intent(in)    ::  a
      real(wp), dimension(nb,nb,*), intent(in)    ::  b
      real(wp), dimension(nb,nb,*), intent(out)   ::  c
      
      integer, dimension(nrow+1),   intent(in)    :: ia
      integer, dimension(*),        intent(in)    :: ja
      integer, dimension(*),        intent(in)    :: ib,jb
      integer, dimension(*),        intent(out)   :: ic,jc

      integer, dimension(*),        intent(inout) :: iw

      integer,                      intent(out)   :: ierr

      integer                  :: k,ka,jj,ii,j,jcol,jpos,kb,len

      real(wp), dimension(nb,nb)                :: scal 

!-----------------------------------------------------------------------
! performs the matrix by matrix product C = A B                         
!-----------------------------------------------------------------------
!
! on entry:  
! ---------                                                             
! nrow  = integer. The row dimension of A = row dimension of C          
! ncol  = integer. The column dimension of B = column dimension of C    
! nb    = integer. The block size of the matrices a, b, and c 
!                                                                       
! a,  ja, ia ; Matrix A in compressed sparse row format.                      
! b,  jb, ib ; Matrix B in compressed sparse row format.                    
!                                                                       
! nzmax = integer. The  length of the arrays c and jc.                  
!         amub will stop if the result matrix C  has a number           
!         of elements that exceeds exceeds nzmax. See ierr.             
!                                                                       
! on return: ! ----------                                                             
! c,  jc, ic ; resulting matrix C in compressed sparse row sparse format.    
!                                                                       
! ierr  = integer. serving as error message.                            
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that amub stopped while computing the       
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! work arrays:                                                          
!------------                                                           
! iw    = integer work array of length equal to the number of columns in A.                                                 

!   The row dimension of B is not needed. However there is no checking  
!   on the condition that ncol(A) = nrow(B).                            

      continue

      len   = 0 
      ic(1) = 1 
      ierr  = 0 

!     initialize array iw.                                              
      do j=1, ncol 
         iw(j) = 0 
      enddo
!                                                                       
      do ii=1, nrow                   !     row i                                                             

         do ka=ia(ii), ia(ii+1)-1 

            scal(:,:) = a (:,:,ka) 
                 jj   = ja(    ka) 

            do kb=ib(jj),ib(jj+1)-1 
               jcol = jb(kb) 
               jpos = iw(jcol) 
               if (jpos  ==  0) then 
                  len = len+1 
                  if (len  >  nzmax) then 
                     ierr = ii 
                     return 
                  endif 
                  jc(len) = jcol 
                  iw(jcol)= len 
                  c(:,:,len)  = matmul( scal(:,:),b(:,:,kb) )
               else 
                  c(:,:,jpos) = c(:,:,jpos) + matmul( scal(:,:),b(:,:,kb) )
               endif 
            enddo

         enddo

         do k=ic(ii), len 
            iw(jc(k)) = 0 
         enddo

         ic(ii+1) = len+1 

      enddo

      end subroutine a_bl_mu_b_bl

!-----------------------------------------------------------------------

  pure subroutine a_bl_mu_b_sc(nrow,ncol,nb,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)                          

      implicit none

      integer,                      intent(in)    :: nrow, ncol, nb, nzmax
      real(wp), dimension(nb,nb,*), intent(in)    ::  a
      real(wp), dimension(*),       intent(in)    ::  b
      real(wp), dimension(nb,nb,*), intent(out)   ::  c
      
      integer, dimension(nrow+1),   intent(in)    :: ia
      integer, dimension(*),        intent(in)    :: ja
      integer, dimension(*),        intent(in)    :: ib,jb
      integer, dimension(*),        intent(out)   :: ic,jc

      integer, dimension(*),        intent(inout) :: iw

      integer,                      intent(out)   :: ierr

      integer                  :: k,ka,jj,ii,j,jcol,jpos,kb,len

      real(wp), dimension(nb,nb)                :: scal 

!-----------------------------------------------------------------------
! performs the matrix by matrix product C = A B                         
!-----------------------------------------------------------------------
!
! on entry:  
! ---------                                                             
! nrow  = integer. The row dimension of A = row dimension of C          
! ncol  = integer. The column dimension of B = column dimension of C    
! nb    = integer. The block size of the matrices a, b, and c 
!                                                                       
! a,  ja, ia ; Matrix A in compressed sparse row format.                      
! b,  jb, ib ; Matrix B in compressed sparse row format.                    
!                                                                       
! nzmax = integer. The  length of the arrays c and jc.                  
!         amub will stop if the result matrix C  has a number           
!         of elements that exceeds exceeds nzmax. See ierr.             
!                                                                       
! on return: ! ----------                                                             
! c,  jc, ic ; resulting matrix C in compressed sparse row sparse format.    
!                                                                       
! ierr  = integer. serving as error message.                            
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that amub stopped while computing the       
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! work arrays:                                                          
!------------                                                           
! iw    = integer work array of length equal to the number of columns in A.                                                 

!   The row dimension of B is not needed. However there is no checking  
!   on the condition that ncol(A) = nrow(B).                            

      continue

      len   = 0 
      ic(1) = 1 
      ierr  = 0 

!     initialize array iw.                                              
      do j=1, ncol 
         iw(j) = 0 
      enddo
!                                                                       
      do ii=1, nrow                   !     row i                                                             

         do ka=ia(ii), ia(ii+1)-1 

            scal(:,:) =  a(:,:,ka) 
                 jj   = ja(    ka) 

            do kb=ib(jj),ib(jj+1)-1 
               jcol = jb(kb) 
               jpos = iw(jcol) 
               if (jpos  ==  0) then 
                  len = len+1 
                  if (len  >  nzmax) then 
                     ierr = ii 
                     return 
                  endif 
                  jc(len) = jcol 
                  iw(jcol)= len 
                  c(:,:,len)  =             + scal(:,:) * b(kb) 
               else 
                  c(:,:,jpos) = c(:,:,jpos) + scal(:,:) * b(kb)
               endif 
            enddo

         enddo

         do k=ic(ii), len 
            iw(jc(k)) = 0 
         enddo

         ic(ii+1) = len+1 

      enddo

      end subroutine a_bl_mu_b_sc

  !-----------------------------------------------------------------------
  subroutine aplb ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
    iw, ierr )

  !*****************************************************************************80
  !
  !! APLB performs the CSR matrix sum C = A + B.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
  !
  !    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
  !
  !    Input, integer ( kind = 4 ) JOB.  When JOB = 0, only the structure
  !    (i.e. the arrays jc, ic) is computed and the
  !    real values are ignored.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! b,
  ! jb,
  ! ib      =  Matrix B in compressed sparse row format.
  !
  ! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number
  !         of elements that exceeds exceeds nzmax. See ierr.
  !
  ! on return:
  !
  ! c,
  ! jc,
  ! ic      = resulting matrix C in compressed sparse row sparse format.
  !
  ! ierr      = integer ( kind = 4 ). serving as error message.
  !         ierr = 0 means normal return,
  !         ierr > 0 means that amub stopped while computing the
  !         i-th row  of C with i = ierr, because the number
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !
  ! iw      = integer ( kind = 4 ) work array of length equal to the number of
  !         columns in A.
  !
    implicit none

    integer :: ncol
    integer :: nrow

    real(wp) a(*)
    real(wp) b(*)
    real(wp) c(*)
    integer :: ia(nrow+1)
    integer :: ib(nrow+1)
    integer :: ic(nrow+1)
    integer :: ierr
    integer :: ii
    integer :: iw(ncol)
    integer :: ja(*)
    integer :: jb(*)
    integer :: jc(*)
    integer :: jcol
    integer :: job
    integer :: jpos
    integer :: k
    integer :: ka
    integer :: kb
    integer :: len
    integer :: nzmax
    logical :: values

    values = ( job /= 0 )
    ierr = 0
    len = 0
    ic(1) = 1
    iw(1:ncol) = 0

    do ii = 1, nrow
  !
  !  Row I.
  !
       do ka = ia(ii), ia(ii+1)-1

          len = len + 1
          jcol = ja(ka)

          if ( nzmax < len ) then
            ierr = ii
            return
          end if

          jc(len) = jcol
          if ( values ) then
            c(len) = a(ka)
          end if
          iw(jcol) = len
       end do

       do kb = ib(ii), ib(ii+1)-1

          jcol = jb(kb)
          jpos = iw(jcol)

          if ( jpos == 0 ) then

             len = len + 1

             if ( nzmax < len ) then
               ierr = ii
               return
             end if

             jc(len) = jcol
             if ( values ) then
               c(len) = b(kb)
             end if
             iw(jcol)= len
          else
             if ( values ) then
               c(jpos) = c(jpos) + b(kb)
             end if
          end if

       end do

       do k = ic(ii), len
         iw(jc(k)) = 0
       end do

       ic(ii+1) = len+1
    end do

    return
  end subroutine aplb

!-----------------------------------------------------------------------

subroutine aplb1 ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, ierr )

!*****************************************************************************80
!
!! APLB1 performs the sum C = A + B for sorted CSR matrices.
!
!  Discussion:
!
!    The difference between this routine and APLB is that here the 
!    resulting matrix is such that the elements of each row are sorted,
!    with increasing column indices in each row, provided the original
!    matrices are sorted in the same way.
!
!    This routine will not work if either of the two input matrices is 
!    not sorted.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
!    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
!    Input, integer ( kind = 4 ) JOB.  When JOB = 0, only the structure

!    (i.e. the arrays jc, ic) is computed and the real values are ignored.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format with entries sorted.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
  implicit none

  integer,                     intent(in)    ::  nrow, ncol, job

  integer,  dimension(nrow+1), intent(in)    :: ia, ib
  integer,  dimension(*),      intent(in)    :: ja, jb
  real(wp), dimension(*),      intent(in)    ::  a,  b

  integer,  dimension(nrow+1), intent(inout) :: ic
  integer,  dimension(*),      intent(inout) :: jc
  real(wp), dimension(*),      intent(inout) ::  c

  integer,                     intent(inout) :: nzmax
  integer,                     intent(out)   :: ierr

  integer :: i
  integer :: j1, j2
  integer :: ka, kamax
  integer :: kb, kbmax
  integer :: kc
  logical :: values
  
  values = ( job /= 0 )
  ierr = 0
  kc = 1
  ic(1) = kc

  do i = 1, nrow

    ka = ia(i)
    kb = ib(i)
    kamax = ia(i+1) - 1
    kbmax = ib(i+1) - 1

    do

      if ( ka <= kamax ) then
        j1 = ja(ka)
      else
        j1 = ncol + 1
      end if

      if ( kb <= kbmax ) then
        j2 = jb(kb)
      else
        j2 = ncol + 1
      end if
!
!  Three cases
!
      if ( j1 == j2 ) then
        if ( values ) then
          c(kc) = a(ka) + b(kb)
        end if
        jc(kc) = j1
        ka = ka + 1
        kb = kb + 1
        kc = kc + 1
      else if ( j1 < j2 ) then
        jc(kc) = j1
        if ( values ) then
          c(kc) = a(ka)
        end if
        ka = ka + 1
        kc = kc + 1
      else if ( j2 < j1 ) then
        jc(kc) = j2
        if ( values ) then
          c(kc) = b(kb)
        end if
        kb = kb + 1
        kc = kc + 1
      end if

      if ( (kc  >  nzmax) .and. (i/=nrow) ) then
        ierr = i
        return
      end if

      if ( kamax < ka .and. kbmax < kb ) then
        exit
      end if

     end do

     ic(i+1) = kc

  end do

  return
end subroutine aplb1

!-----------------------------------------------------------------------

subroutine a_bl_pl_b_bl1( nrow, ncol, nblk, a, ja, ia, b, jb, ib,  &
                                          & c, jc, ic, nzmax, ierr )

!*****************************************************************************80
!
!! APLB1 performs the sum C = A + B for sorted CSR matrices.
!
!  Discussion:
!
!    The difference between this routine and APLB is that here the 
!    resulting matrix is such that the elements of each row are sorted,
!    with increasing column indices in each row, provided the original
!    matrices are sorted in the same way.
!
!    This routine will not work if either of the two input matrices is 
!    not sorted.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of A and B.
!    Input, integer ( kind = 4 ) NCOL, the column dimension of A and B.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format with entries sorted.
!
!  b,
! jb,
! ib      =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax      = integer ( kind = 4 ). The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
!  c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr      = integer ( kind = 4 ). serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
  implicit none

  integer,                          intent(in)    ::  nrow, ncol, nblk, nzmax

  integer,  dimension(nrow+1),      intent(in)    :: ia, ib
  integer,  dimension(*),           intent(in)    :: ja, jb
  real(wp), dimension(nblk,nblk,*), intent(in)    ::  a,  b

  integer,  dimension(nrow+1),      intent(inout) :: ic
  integer,  dimension(*),           intent(inout) :: jc
  real(wp), dimension(nblk,nblk,*), intent(inout) ::  c

  integer,                          intent(out)   :: ierr

  integer :: i
  integer :: j1, j2
  integer :: ka, kamax
  integer :: kb, kbmax
  integer :: kc

  ierr = 0
  kc = 1
  ic(1) = kc

  do i = 1, nrow

       ka = ia(i)
       kb = ib(i)
    kamax = ia(i+1) - 1
    kbmax = ib(i+1) - 1

    do

      if ( ka <= kamax ) then
        j1 = ja(ka)
      else
        j1 = ncol + 1
      end if

      if ( kb <= kbmax ) then
        j2 = jb(kb)
      else
        j2 = ncol + 1
      end if
!
!  Three cases
!
      if ( j1 == j2 ) then
        c(:,:,kc) = a(:,:,ka) + b(:,:,kb)
        jc(kc) = j1
        ka = ka + 1
        kb = kb + 1
        kc = kc + 1
      else if ( j1 < j2 ) then
        jc(kc) = j1
        c(:,:,kc) = a(:,:,ka)
        ka = ka + 1
        kc = kc + 1
      else if ( j2 < j1 ) then
        jc(kc) = j2
        c(:,:,kc) = b(:,:,kb)
        kb = kb + 1
        kc = kc + 1
      end if

      if ( (kc  >  nzmax) .and. (i/=nrow) ) then
        ierr = i
        return
      end if

      if ( kamax < ka .and. kbmax < kb ) then
        exit
      end if

     end do

     ic(i+1) = kc

  end do

  return
end subroutine a_bl_pl_b_bl1

!-----------------------------------------------------------------------

      subroutine aplsb (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,ierr)
      integer nrow,ncol,nzmax,ierr
      real(wp) a(*), b(*), c(*), s 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1) 
!-----------------------------------------------------------------------
! performs the operation C = A+s B for matrices in sorted CSR format.   
! the difference with aplsb is that the resulting matrix is such that   
! the elements of each row are sorted with increasing column indices in 
! each row, provided the original matrices are sorted in the same way.  
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A and B                   
! ncol  = integer. The column dimension of A and B.                     
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format with entries sorted   
!                                                                       
! s        = real. scalar factor for B.                                 
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        =  Matrix B in compressed sparse row format with entries sor
!        ascendly in each row                                           
!                                                                       
! nzmax        = integer. The  length of the arrays c and jc.           
!         amub will stop if the result matrix C  has a number           
!         of elements that exceeds exceeds nzmax. See ierr.             
!                                                                       
! on return:                                                            
!----------                                                             
! c,                                                                    
! jc,                                                                   
! ic        = resulting matrix C in compressed sparse row sparse format 
!         with entries sorted ascendly in each row.                     
!                                                                       
! ierr        = integer. serving as error message.                      
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that amub stopped while computing the       
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! Notes:                                                                
!-------                                                                
!     this will not work if any of the two input matrices is not sorted 
!-----------------------------------------------------------------------
      integer i, j1, j2, ka, kamax, kb, kbmax, kc

      ierr = 0 
      kc = 1 
      ic(1) = kc 
!                                                                       
!     the following loop does a merge of two sparse rows + adds  them.  
!                                                                       
      do 6 i=1, nrow 
         ka = ia(i) 
         kb = ib(i) 
         kamax = ia(i+1)-1 
         kbmax = ib(i+1)-1 
    5    continue 
!                                                                       
!     this is a while  -- do loop --                                    
!                                                                       
         if (ka  <=  kamax .or. kb  <=  kbmax) then 
!                                                                       
            if (ka  <=  kamax) then 
               j1 = ja(ka) 
            else 
!     take j1 large enough  that always j2  <  j1                      
               j1 = ncol+1 
            endif 
            if (kb  <=  kbmax) then 
               j2 = jb(kb) 
            else 
!     similarly take j2 large enough  that always j1  <  j2            
               j2 = ncol+1 
            endif 
!                                                                       
!     three cases                                                       
!                                                                       
            if (j1  ==  j2) then 
               c(kc) = a(ka)+s*b(kb) 
               jc(kc) = j1 
               ka = ka+1 
               kb = kb+1 
               kc = kc+1 
            else if (j1  <  j2) then 
               jc(kc) = j1 
               c(kc) = a(ka) 
               ka = ka+1 
               kc = kc+1 
            else if (j1  >  j2) then 
               jc(kc) = j2 
               c(kc) = s*b(kb) 
               kb = kb+1 
               kc = kc+1 
            endif 
            if (kc  >  nzmax) goto 999 
            goto 5 
!                                                                       
!     end while loop                                                    
!                                                                       
         endif 
         ic(i+1) = kc 
    6 continue 
      return 
  999 ierr = i 
      return 
!------------end-of-aplsb --------------------------------------------- 
!-----------------------------------------------------------------------
      end subroutine aplsb

!-----------------------------------------------------------------------
      subroutine aplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,           &
     &     nzmax,ierr)                                                  
     integer nrow, ncol, nzmax, ierr 
     real(wp) a(*), b(*), c(*), s 
     integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1) 
!-----------------------------------------------------------------------
! performs the operation C = A+s B for matrices in sorted CSR format.   
! the difference with aplsb is that the resulting matrix is such that   
! the elements of each row are sorted with increasing column indices in 
! each row, provided the original matrices are sorted in the same way.  
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A and B                   
! ncol  = integer. The column dimension of A and B.                     
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format with entries sorted   
!                                                                       
! s        = real. scalar factor for B.                                 
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        =  Matrix B in compressed sparse row format with entries sor
!        ascendly in each row                                           
!                                                                       
! nzmax        = integer. The  length of the arrays c and jc.           
!         amub will stop if the result matrix C  has a number           
!         of elements that exceeds exceeds nzmax. See ierr.             
!                                                                       
! on return:                                                            
!----------                                                             
! c,                                                                    
! jc,                                                                   
! ic        = resulting matrix C in compressed sparse row sparse format 
!         with entries sorted ascendly in each row.                     
!                                                                       
! ierr        = integer. serving as error message.                      
!         ierr = 0 means normal return,                                 
!         ierr  >  0 means that amub stopped while computing the       
!         i-th row  of C with i=ierr, because the number                
!         of elements in C exceeds nzmax.                               
!                                                                       
! Notes:                                                                
!-------                                                                
!     this will not work if any of the two input matrices is not sorted 
!-----------------------------------------------------------------------
      integer i, j1, j2, ka, kamax, kb, kbmax, kc

      ierr = 0 
      kc = 1 
      ic(1) = kc 
!                                                                       
!     the following loop does a merge of two sparse rows + adds  them.  
!                                                                       
      do 6 i=1, nrow 
         ka = ia(i) 
         kb = ib(i) 
         kamax = ia(i+1)-1 
         kbmax = ib(i+1)-1 
    5    continue 
!                                                                       
!     this is a while  -- do loop --                                    
!                                                                       
         if (ka  <=  kamax .or. kb  <=  kbmax) then 
!                                                                       
            if (ka  <=  kamax) then 
               j1 = ja(ka) 
            else 
!     take j1 large enough  that always j2  <  j1                      
               j1 = ncol+1 
            endif 
            if (kb  <=  kbmax) then 
               j2 = jb(kb) 
            else 
!     similarly take j2 large enough  that always j1  <  j2            
               j2 = ncol+1 
            endif 
!                                                                       
!     three cases                                                       
!                                                                       
            if (j1  ==  j2) then 
               c(kc) = a(ka)+s*b(kb) 
               jc(kc) = j1 
               ka = ka+1 
               kb = kb+1 
               kc = kc+1 
            else if (j1  <  j2) then 
               jc(kc) = j1 
               c(kc) = a(ka) 
               ka = ka+1 
               kc = kc+1 
            else if (j1  >  j2) then 
               jc(kc) = j2 
               c(kc) = s*b(kb) 
               kb = kb+1 
               kc = kc+1 
            endif 
            if (kc  >  nzmax) goto 999 
            goto 5 
!                                                                       
!     end while loop                                                    
!                                                                       
         endif 
         ic(i+1) = kc 
    6 continue 
      return 
  999 ierr = i 
      return 
!------------end-of-aplsb1 ---------------------------------------------
!-----------------------------------------------------------------------
      end subroutine aplsb1

!-----------------------------------------------------------------------
!      subroutine apmbt (nrow,ncol,job,a,ja,ia,b,jb,ib,                  &
!     &     c,jc,ic,nzmax,iw,ierr)                                       
!      real(wp) a(*), b(*), c(*) 
!      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),iw(*) 
!!-----------------------------------------------------------------------
!! performs the matrix sum  C = A + transp(B) or C = A - transp(B)       
!!-----------------------------------------------------------------------
!! on entry:                                                             
!! ---------                                                             
!! nrow        = integer. The row dimension of A and transp(B)           
!! ncol  = integer. The column dimension of A. Also the row              
!!                  dimension of B.                                      
!!                                                                       
!! job        = integer. if job = -1, apmbt will compute C= A - transp(B)
!!         (structure + values)                                          
!!         if (job  ==  1)  it will compute C=A+transp(A)                
!!         (structure+ values)                                           
!!         if (job  ==  0) it will compute the structure of              
!!         C= A+/-transp(B) only (ignoring all real values).             
!!         any other value of job will be treated as  job=1              
!! a,                                                                    
!! ja,                                                                   
!! ia    = Matrix A in compressed sparse row format.                     
!!                                                                       
!! b,                                                                    
!! jb,                                                                   
!! ib        =  Matrix B in compressed sparse row format.                
!!                                                                       
!! nzmax        = integer. The  length of the arrays c, jc, and ic.      
!!         amub will stop if the result matrix C  has a number           
!!         of elements that exceeds exceeds nzmax. See ierr.             
!!                                                                       
!! on return:                                                            
!!----------                                                             
!! c,                                                                    
!! jc,                                                                   
!! ic        = resulting matrix C in compressed sparse row format.       
!!                                                                       
!! ierr        = integer. serving as error message.                      
!!         ierr = 0 means normal return.                                 
!!         ierr = -1 means that nzmax was  <  either the number of      
!!         nonzero elements of A or the number of nonzero elements in B. 
!!         ierr  >  0 means that amub stopped while computing the       
!!         i-th row  of C with i=ierr, because the number                
!!         of elements in C exceeds nzmax.                               
!!                                                                       
!! work arrays:                                                          
!!------------                                                           
!! iw        = integer work array of length at least max(ncol,nrow)      
!!                                                                       
!! Notes:                                                                
!!------- It is important to note that here all of three arrays c, ic,   
!!        and jc are assumed to be of length nnz(c). This is because     
!!        the matrix is internally converted in coordinate format.       
!!                                                                       
!!-----------------------------------------------------------------------
!      logical values 
!      values = (job  /=  0) 
!!                                                                       
!      ierr = 0 
!      do 1 j=1, ncol 
!         iw(j) = 0 
!    1 continue 
!!                                                                       
!      nnza = ia(nrow+1)-1 
!      nnzb = ib(ncol+1)-1 
!      len = nnzb 
!      if (nzmax  <  nnzb .or. nzmax  <  nnza) then 
!         ierr = -1 
!         return 
!      endif 
!!                                                                       
!! trasnpose matrix b into c                                             
!!                                                                       
!      ljob = 0 
!      if (values) ljob = 1 
!      ipos = 1 
!      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
!!-----------------------------------------------------------------------
!      if (job  ==  -1) then 
!         do 2 k=1,len 
!            c(k) = -c(k) 
!    2    continue 
!      endif 
!!                                                                       
!!--------------- main loop -------------------------------------------- 
!!                                                                       
!      do 500 ii=1, nrow 
!         do 200 k = ic(ii),ic(ii+1)-1 
!            iw(jc(k)) = k 
!  200    continue 
!!-----------------------------------------------------------------------
!         do 300 ka = ia(ii), ia(ii+1)-1 
!            jcol = ja(ka) 
!            jpos = iw(jcol) 
!            if (jpos  ==  0) then 
!!                                                                       
!!     if fill-in append in coordinate format to matrix.                 
!!                                                                       
!               len = len+1 
!               if (len  >  nzmax) goto 999 
!               jc(len) = jcol 
!                                                                        
!               ic(len) = ii 
!               if (values) c(len)  = a(ka) 
!            else 
!!     else do addition.                                                 
!               if (values) c(jpos) = c(jpos) + a(ka) 
!            endif 
!  300    continue 
!         do 301 k=ic(ii), ic(ii+1)-1 
!            iw(jc(k)) = 0 
!  301    continue 
!  500 continue 
!!                                                                       
!!     convert first part of matrix (without fill-ins) into coo format   
!!                                                                       
!      ljob = 2 
!      if (values) ljob = 3 
!      do 501 i=1, nrow+1 
!         iw(i) = ic(i) 
!  501 continue 
!      call csrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr) 
!!                                                                       
!!     convert the whole thing back to csr format.                       
!!                                                                       
!      ljob = 0 
!      if (values) ljob = 1 
!      call coicsr (nrow,len,ljob,c,jc,ic,iw) 
!      return 
!  999 ierr = ii 
!      return 
!!--------end-of-apmbt---------------------------------------------------
!!-----------------------------------------------------------------------
!      end subroutine apmbt
!
!-----------------------------------------------------------------------
!      subroutine aplsbt(nrow,ncol,a,ja,ia,s,b,jb,ib,                    &
!     &     c,jc,ic,nzmax,iw,ierr)                                       
!      real(wp) a(*), b(*), c(*), s 
!      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),iw(*) 
!!-----------------------------------------------------------------------
!! performs the matrix sum  C = A + transp(B).                           
!!-----------------------------------------------------------------------
!! on entry:                                                             
!! ---------                                                             
!! nrow        = integer. The row dimension of A and transp(B)           
!! ncol  = integer. The column dimension of A. Also the row              
!!                  dimension of B.                                      
!!                                                                       
!! a,                                                                    
!! ja,                                                                   
!! ia    = Matrix A in compressed sparse row format.                     
!!                                                                       
!! s        = real. scalar factor for B.                                 
!!                                                                       
!!                                                                       
!! b,                                                                    
!! jb,                                                                   
!! ib        =  Matrix B in compressed sparse row format.                
!!                                                                       
!! nzmax        = integer. The  length of the arrays c, jc, and ic.      
!!         amub will stop if the result matrix C  has a number           
!!         of elements that exceeds exceeds nzmax. See ierr.             
!!                                                                       
!! on return:                                                            
!!----------                                                             
!! c,                                                                    
!! jc,                                                                   
!! ic        = resulting matrix C in compressed sparse row format.       
!!                                                                       
!! ierr        = integer. serving as error message.                      
!!         ierr = 0 means normal return.                                 
!!         ierr = -1 means that nzmax was  <  either the number of      
!!         nonzero elements of A or the number of nonzero elements in B. 
!!         ierr  >  0 means that amub stopped while computing the       
!!         i-th row  of C with i=ierr, because the number                
!!         of elements in C exceeds nzmax.                               
!!                                                                       
!! work arrays:                                                          
!!------------                                                           
!! iw        = integer work array of length at least max(nrow,ncol)      
!!                                                                       
!! Notes:                                                                
!!------- It is important to note that here all of three arrays c, ic,   
!!        and jc are assumed to be of length nnz(c). This is because     
!!        the matrix is internally converted in coordinate format.       
!!                                                                       
!!-----------------------------------------------------------------------
!      ierr = 0 
!      do 1 j=1, ncol 
!         iw(j) = 0 
!    1 continue 
!!                                                                       
!      nnza = ia(nrow+1)-1 
!      nnzb = ib(ncol+1)-1 
!      len = nnzb 
!      if (nzmax  <  nnzb .or. nzmax  <  nnza) then 
!         ierr = -1 
!         return 
!      endif 
!!                                                                       
!!     transpose matrix b into c                                         
!!                                                                       
!      ljob = 1 
!      ipos = 1 
!      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
!      do 2 k=1,len 
!    2    c(k) = c(k)*s 
!!                                                                       
!!     main loop. add rows from ii = 1 to nrow.                          
!!                                                                       
!         do 500 ii=1, nrow 
!!     iw is used as a system to recognize whether there                 
!!     was a nonzero element in c.                                       
!            do 200 k = ic(ii),ic(ii+1)-1 
!               iw(jc(k)) = k 
!  200       continue 
!!                                                                       
!            do 300 ka = ia(ii), ia(ii+1)-1 
!               jcol = ja(ka) 
!               jpos = iw(jcol) 
!           if (jpos  ==  0) then 
!!                                                                       
!!     if fill-in append in coordinate format to matrix.                 
!!                                                                       
!              len = len+1 
!              if (len  >  nzmax) goto 999 
!              jc(len) = jcol 
!              ic(len) = ii 
!              c(len)  = a(ka) 
!           else 
!!     else do addition.                                                 
!              c(jpos) = c(jpos) + a(ka) 
!           endif 
!  300   continue 
!        do 301 k=ic(ii), ic(ii+1)-1 
!           iw(jc(k)) = 0 
!  301   continue 
!  500 continue 
!!                                                                       
!!     convert first part of matrix (without fill-ins) into coo format   
!!                                                                       
!      ljob = 3 
!      do 501 i=1, nrow+1 
!         iw(i) = ic(i) 
!  501 continue 
!      call csrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr) 
!!                                                                       
!!     convert the whole thing back to csr format.                       
!!                                                                       
!      ljob = 1 
!      call coicsr (nrow,len,ljob,c,jc,ic,iw) 
!      return 
!  999 ierr = ii 
!      return 
!!--------end-of-aplsbt--------------------------------------------------
!!-----------------------------------------------------------------------
!      end subroutine aplsbt

!-----------------------------------------------------------------------
      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib) 
      integer nrow, job 
      real(wp) a(*), b(*), diag(nrow), scal 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
!-----------------------------------------------------------------------
! performs the matrix by matrix product B = Diag * A  (in place)        
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! job   = integer. job indicator. Job=0 means get array b only          
!         job = 1 means get b, and the integer arrays ib, jb.           
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! diag = diagonal matrix stored as a vector dig(1:n)                    
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        = resulting matrix B in compressed sparse row sparse format.
!                                                                       
! Notes:                                                                
!-------                                                                
! 1)        The column dimension of A is not needed.                    
! 2)        algorithm in place (B can take the place of A).             
!           in this case use job=0.                                     
!-----------------------------------------------------------------      
      integer ii, k, k1, k2
      do 1 ii=1,nrow 
!                                                                       
!     normalize each row                                                
!                                                                       
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         scal = diag(ii) 
         do 2 k=k1, k2 
            b(k) = a(k)*scal 
    2    continue 
    1 continue 
!                                                                       
      if (job  ==  0) return 
!                                                                       
      do 3 ii=1, nrow+1 
         ib(ii) = ia(ii) 
    3 continue 
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k) 
   31 continue 

      return 
      end subroutine diamua

!-----------------------------------------------------------------------
      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib) 
      integer nrow, job
      real(wp) a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)
      

!-----------------------------------------------------------------------
! performs the matrix by matrix product B = A * Diag  (in place)        
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! job   = integer. job indicator. Job=0 means get array b only          
!         job = 1 means get b, and the integer arrays ib, jb.           
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! diag = diagonal matrix stored as a vector dig(1:n)                    
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        = resulting matrix B in compressed sparse row sparse format.
!                                                                       
! Notes:                                                                
!-------                                                                
! 1)        The column dimension of A is not needed.                    
! 2)        algorithm in place (B can take the place of A).             
!-----------------------------------------------------------------      
      integer ii, k, k1, k2

      do 1 ii=1,nrow 
!                                                                       
!     scale each element                                                
!                                                                       
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         do 2 k=k1, k2 
            b(k) = a(k)*diag(ja(k)) 
    2    continue 
    1 continue 
!                                                                       
      if (job  ==  0) return 
!                                                                       
      do 3 ii=1, nrow+1 
         ib(ii) = ia(ii) 
    3 continue 
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k) 
   31 continue 
      return 

      end subroutine amudia

!-----------------------------------------------------------------------

      subroutine aplsca (nrow, a, ja, ia, scal,iw) 
      integer nrow
      real(wp) a(*), scal 
      integer ja(*), ia(nrow+1),iw(*) 
!-----------------------------------------------------------------------
! Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I  
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia    = Matrix A in compressed sparse row format.                     
!                                                                       
! scal  = real. scalar to add to the diagonal entries.                  
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia        = matrix A with diagonal elements shifted (or created).     
!                                                                       
! iw    = integer work array of length n. On return iw will             
!         contain  the positions of the diagonal entries in the         
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,          
!         are the values/column indices of the diagonal elements        
!         of the output matrix. ).                                      
!                                                                       
! Notes:                                                                
!-------                                                                
!     The column dimension of A is not needed.                          
!     important: the matrix a may be expanded slightly to allow for     
!     additions of nonzero elements to previously nonexisting diagonals.
!     The is no checking as to whether there is enough space appended   
!     to the arrays a and ja. if not sure allow for n additional        
!     elemnts.                                                          
!     coded by Y. Saad. Latest version July, 19, 1990                   
!-----------------------------------------------------------------------
      logical test
      integer icount, ii, j, k, k1, k2, ko
!                                                                       
      call diapos (nrow,ja,ia,iw) 
      icount = 0 
      do 1 j=1, nrow 
         if (iw(j)  ==  0) then 
            icount = icount+1 
         else 
            a(iw(j)) = a(iw(j)) + scal 
         endif 
    1 continue 
!                                                                       
!     if no diagonal elements to insert in data structure return.       
!                                                                       
      if (icount  ==  0) return 
!                                                                       
! shift the nonzero elements if needed, to allow for created            
! diagonal elements.                                                    
!                                                                       
      ko = ia(nrow+1)+icount 
!                                                                       
!     copy rows backward                                                
!                                                                       
      do 5 ii=nrow, 1, -1 
!                                                                       
!     go through  row ii                                                
!                                                                       
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         ia(ii+1) = ko 
         test = (iw(ii)  ==  0) 
         do 4 k = k2,k1,-1 
            j = ja(k) 
            if (test .and. (j  <  ii)) then 
               test = .false. 
               ko = ko - 1 
               a(ko) = scal 
               ja(ko) = ii 
               iw(ii) = ko 
            endif 
            ko = ko-1 
            a(ko) = a(k) 
            ja(ko) = j 
    4    continue 
!     diagonal element has not been added yet.                          
         if (test) then 
            ko = ko-1 
            a(ko) = scal 
            ja(ko) = ii 
            iw(ii) = ko 
         endif 
    5 continue 
      ia(1) = ko 

      return 
      end subroutine aplsca

!-----------------------------------------------------------------------

      subroutine a_bl_pl_b_bl_scal(nrow, nb, a, ja, ia, scal, iw) 

      implicit none

      integer,                         intent(in)    :: nrow, nb
      real(wp),                        intent(in)    :: scal

      real(wp), dimension(nb,nb,*),    intent(inout) ::  a

      integer, dimension(nrow+1),      intent(inout) :: ia
      integer, dimension(*),           intent(inout) :: ja

      integer, dimension(*),           intent(inout) :: iw

      real(wp), dimension(nb,nb)                     :: eye

      integer          :: i,icount, ii,j,k,k1,k2,ko
!-----------------------------------------------------------------------
! Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I  
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow        = integer. The row dimension of A                         
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia    = Matrix A in compressed sparse row format.                     
!                                                                       
! scal  = real. scalar to add to the diagonal entries.                  
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia        = matrix A with diagonal elements shifted (or created).     
!                                                                       
! iw    = integer work array of length n. On return iw will             
!         contain  the positions of the diagonal entries in the         
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,          
!         are the values/column indices of the diagonal elements        
!         of the output matrix. ).                                      
!                                                                       
! Notes:                                                                
!-------                                                                
!     The column dimension of A is not needed.                          
!     important: the matrix a may be expanded slightly to allow for     
!     additions of nonzero elements to previously nonexisting diagonals.
!     The is no checking as to whether there is enough space appended   
!     to the arrays a and ja. if not sure allow for n additional        
!     elemnts.                                                          
!     coded by Y. Saad. Latest version July, 19, 1990                   
!-----------------------------------------------------------------------
      logical test 

      eye = 0.0_wp
      do i = 1,nb
        eye(i,i) = 1.0_wp
      enddo
!                                                                       
      call diapos (nrow,ja,ia,iw) 
      icount = 0 

      do 1 j=1, nrow 
         if (iw(j)  ==  0) then 
            icount = icount+1 
         else 
            a(:,:,iw(j)) = a(:,:,iw(j)) + scal * eye(:,:)
         endif 
    1 continue 
!                                                                       
!     if no diagonal elements to insert in data structure return.       
!                                                                       
      if (icount  ==  0) return 
!                                                                       
! shift the nonzero elements if needed, to allow for created            
! diagonal elements.                                                    
!                                                                       
      ko = ia(nrow+1)+icount 
!                                                                       
!     copy rows backward                                                
!                                                                       
      do 5 ii=nrow, 1, -1 
!                                                                       
!     go through  row ii                                                
!                                                                       
         k1 = ia(ii) 
         k2 = ia(ii+1)-1 
         ia(ii+1) = ko 
         test = (iw(ii)  ==  0) 
         do 4 k = k2,k1,-1 
            j = ja(k) 
            if (test .and. (j  <  ii)) then 
               test = .false. 
               ko = ko - 1 
               a(:,:,ko) = scal * eye(:,:)
               ja(ko) = ii 
               iw(ii) = ko 
            endif 
            ko = ko-1 
            a(:,:,ko) = a(:,:,k) 
            ja(ko) = j 
    4    continue 
!     diagonal element has not been added yet.                          
         if (test) then 
            ko = ko-1 
            a(:,:,ko) = scal * eye(:,:)
            ja(ko) = ii 
            iw(ii) = ko 
         endif 
    5 continue 
      ia(1) = ko 

      return 
      end subroutine a_bl_pl_b_bl_scal

!-----------------------------------------------------------------------

      subroutine apldia (nrow, job, a, ja, ia, diag, b, jb, ib, iw) 
      integer nrow
      real(wp) a(*), b(*), diag(nrow) 
      integer job,ja(*),jb(*), ia(nrow+1),ib(nrow+1), iw(*) 
!-----------------------------------------------------------------------
! Adds a diagonal matrix to a general sparse matrix:  B = A + Diag      
!-----------------------------------------------------------------------
! on entry:                                                             
! ---------                                                             
! nrow  = integer. The row dimension of A                         
!                                                                       
! job   = integer. job indicator. Job=0 means get array b only          
!         (i.e. assume that a has already been copied into array b,     
!         or that algorithm is used in place. ) For all practical       
!         purposes enter job=0 for an in-place call and job=1 otherwise 
!                                                                       
!         Note: in case there are missing diagonal elements in A,       
!         then the option job =0 will be ignored, since the algorithm   
!         must modify the data structure (i.e. jb, ib) in this          
!         situation.                                                    
!                                                                       
! a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! diag = diagonal matrix stored as a vector dig(1:n)                    
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
! b,                                                                    
! jb,                                                                   
! ib        = resulting matrix B in compressed sparse row sparse format.
!                                                                       
!                                                                       
! iw    = integer work array of length n. On return iw will             
!         contain  the positions of the diagonal entries in the         
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,          
!         are the values/column indices of the diagonal elements        
!         of the output matrix. ).                                      
!                                                                       
! Notes:                                                                
!-------                                                                
! 1)        The column dimension of A is not needed.                    
! 2)        algorithm in place (b, jb, ib, can be the same as           
!           a, ja, ia, on entry). See comments for parameter job.       
!                                                                       
! coded by Y. Saad. Latest version July, 19, 1990                       
!-----------------------------------------------------------------      
      logical test
      integer icount, ii, j, k, k1, k2, ko, nnz
!                                                                       
!     copy integer arrays into b's data structure if required           
!                                                                       
      if (job  /=  0) then 
         nnz = ia(nrow+1)-1 
         do 2  k=1, nnz 
            jb(k) = ja(k) 
            b(k)  = a(k) 
    2    continue 
         do 3 k=1, nrow+1 
            ib(k) = ia(k) 
    3    continue 
      endif 
                                                                        
!     get positions of diagonal elements in data structure.             
                                                                        
      call diapos (nrow,ja,ia,iw) 
                                                                        
!     count number of holes in diagonal and add diag(*) elements to     
!     valid diagonal entries.                                           
                                                                        
      icount = 0 
      do 1 j=1, nrow 
         if (iw(j)  ==  0) then 
            icount = icount+1 
         else 
            b(iw(j)) = a(iw(j)) + diag(j) 
         endif 
    1 continue 
                                                                       
!     if no diagonal elements to insert return                          
                                                                        
      if (icount  ==  0) return 
!                                                                       
!     shift the nonzero elements if needed, to allow for created        
!     diagonal elements.                                                
!                                                                       
      ko = ib(nrow+1)+icount 
                                                                        
!     copy rows backward                                                
      do 5 ii=nrow, 1, -1 
!        go through  row ii                                                
         k1 = ib(ii) 
         k2 = ib(ii+1)-1 
         ib(ii+1) = ko 
         test = (iw(ii)  ==  0) 
         do 4 k = k2,k1,-1 
            j = jb(k) 
            if (test .and. (j  <  ii)) then 
               test = .false. 
               ko = ko - 1 
               b(ko) = diag(ii) 
               jb(ko) = ii 
               iw(ii) = ko 
            endif 
            ko = ko-1 
            b(ko) = a(k) 
            jb(ko) = j 
    4    continue 
!     diagonal element has not been added yet.                          
         if (test) then 
            ko = ko-1 
            b(ko) =  diag(ii) 
            jb(ko) = ii 
            iw(ii) = ko 
         endif 
    5 continue 
      ib(1) = ko 

      return 
      end subroutine apldia

!=============================================================================80

      subroutine a_bl_pl_b_bl_diag(nrow, nb, job, a, ja, ia, diag, b, jb, ib, iw) 

      implicit none

      integer,                         intent(in)    :: nrow, nb, job
      real(wp), dimension(nb,nb,*),    intent(in)    ::  a
      real(wp), dimension(nb,nb,nrow), intent(in)    ::  diag

      real(wp), dimension(nb,nb,*),    intent(inout) ::  b

      integer, dimension(nrow+1),      intent(in)    :: ia
      integer, dimension(*),           intent(in)    :: ja

      integer, dimension(nrow+1),      intent(inout) :: ib
      integer, dimension(*),           intent(inout) :: jb

      integer, dimension(*),           intent(inout) :: iw

      integer              :: ii,j,k, icount, k1,k2,ko, nnz
      logical test 

!-----------------------------------------------------------------------
! Adds a diagonal matrix to a general sparse matrix:  B = A + Diag      
!-----------------------------------------------------------------------

! on entry:                                                             
! ---------                                                             
! nrow  = integer. The row dimension of A                         
!                                                                       
! job   = integer. job indicator. Job=0 means get array b only          
!         (i.e. assume that a has already been copied into array b,     
!         or that algorithm is used in place. ) For all practical       
!         purposes enter job=0 for an in-place call and job=1 otherwise 
!                                                                       
!         Note: in case there are missing diagonal elements in A,       
!         then the option job =0 will be ignored, since the algorithm   
!         must modify the data structure (i.e. jb, ib) in this          
!         situation.                                                    
!                                                                       
!  a,                                                                    
! ja,                                                                   
! ia   = Matrix A in compressed sparse row format.                      
!                                                                       
! diag = diagonal matrix stored as a vector dig(1:n)                    
!                                                                       
! on return:                                                            
!----------                                                             
!                                                                       
!  b,                                                                    
! jb,                                                                   
! ib        = resulting matrix B in compressed sparse row sparse format.
!                                                                       
!                                                                       
! iw    = integer work array of length n. On return iw will             
!         contain  the positions of the diagonal entries in the         
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,          
!         are the values/column indices of the diagonal elements        
!         of the output matrix. ).                                      
!                                                                       
! Notes:                                                                
!-------                                                                
! 1)        The column dimension of A is not needed.                    
! 2)        algorithm in place (b, jb, ib, can be the same as           
!           a, ja, ia, on entry). See comments for parameter job.       
!                                                                       
! coded by Y. Saad. Latest version July, 19, 1990                       
!-----------------------------------------------------------------      
!                                                                       
!     copy integer arrays into b's data structure if required           
!                                                                       
      if (job  /=  0) then 
         nnz = ia(nrow+1)-1 
         do 2  k=1, nnz 
            jb(k)      = ja(k) 
             b(:,:,k)  = a(:,:,k) 
    2    continue 
         do 3 k=1, nrow+1 
            ib(k) = ia(k) 
    3    continue 
      endif 
                                                                        
!     get positions of diagonal elements in data structure.             
                                                                        
      call diapos (nrow,ja,ia,iw) 
                                                                        
!     count number of holes in diagonal and add diag(*) elements to     
!     valid diagonal entries.                                           
                                                                        
      icount = 0 
      do 1 j=1, nrow 
         if (iw(j)  ==  0) then 
            icount = icount+1 
         else 
            b(:,:,iw(j)) = a(:,:,iw(j)) + diag(:,:,j) 
         endif 
    1 continue 
                                                                       
!     if no diagonal elements to insert return                          
                                                                        
      if (icount  ==  0) return 
!                                                                       
!     shift the nonzero elements if needed, to allow for created        
!     diagonal elements.                                                
!                                                                       
      ko = ib(nrow+1)+icount 
                                                                        
!     copy rows backward                                                
      do 5 ii=nrow, 1, -1 
!        go through  row ii                                                
         k1 = ib(ii) 
         k2 = ib(ii+1)-1 
         ib(ii+1) = ko 
         test = (iw(ii)  ==  0) 
         do 4 k = k2,k1,-1 
            j = jb(k) 
            if (test .and. (j  <  ii)) then 
               test = .false. 
               ko = ko - 1 
               b(:,:,ko) = diag(:,:,ii) 
               jb(ko) = ii 
               iw(ii) = ko 
            endif 
            ko = ko-1 
            b(:,:,ko) = a(:,:,k) 
            jb(ko) = j 
    4    continue 
!     diagonal element has not been added yet.                          
         if (test) then 
            ko = ko-1 
            b(:,:,ko) =  diag(:,:,ii) 
            jb(ko) = ii 
            iw(ii) = ko 
         endif 
    5 continue 
      ib(1) = ko 

      return 
      end subroutine a_bl_pl_b_bl_diag

!=============================================================================80
 
      subroutine qsortd(x,ind,n)

!       Code converted using TO_F90 by Alan Miller
!       Date: 2002-12-18  Time: 11:55:47

      real(wp), INTENT(IN)    :: x(:)
      INTEGER,  INTENT(OUT)   :: ind(:)
      INTEGER,  INTENT(IN)    :: n

!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL(wp)
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

      INTEGER     :: iu(21), il(21)
      INTEGER     :: m, i, j, k, l, ij, it, itt, indx
      REAL        :: r
      REAL(wp)    :: t

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

      IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

      DO  i = 1, n
        ind(i) = i
      END DO
      m = 1
      i = 1
      j = n
      r = .375

! TOP OF LOOP

      20 IF (i >= j) GO TO 70
      IF (r <= .5898437) THEN
        r = r + .0390625
      ELSE
        r = r - .21875
      END IF

! INITIALIZE K

      30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

      ij = i + r*(j-i)
      it = ind(ij)
      t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) > t) THEN
        ind(ij) = indx
        ind(i) = it
        it = indx
        t = x(it)
      END IF

! INITIALIZE L

      l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

      indx = ind(j)
      IF (x(indx) >= t) GO TO 50
      ind(ij) = indx
      ind(j) = it
      it = indx
      t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) <= t) GO TO 50
      ind(ij) = indx
      ind(i) = it
      it = indx
      t = x(it)
      GO TO 50

! INTERCHANGE ELEMENTS K AND L

      40 itt = ind(l)
      ind(l) = ind(k)
      ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

      50 l = l - 1
      indx = ind(l)
      IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

      60 k = k + 1
      indx = ind(k)
      IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

      IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

      IF (l-i > j-k) THEN
        il(m) = i
        iu(m) = l
        i = k
        m = m + 1
        GO TO 80
      END IF

      il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

      70 m = m - 1
      IF (m == 0) RETURN
      i = il(m)
      j = iu(m)

      80 IF (j-i >= 11) GO TO 30
      IF (i == 1) GO TO 20
      i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

      90 i = i + 1
      IF (i == j) GO TO 70
      indx = ind(i+1)
      t = x(indx)
      it = indx
      indx = ind(i)
      IF (x(indx) <= t) GO TO 90
      k = i

      100 ind(k+1) = ind(k)
      k = k - 1
      indx = ind(k)
      IF (t < x(indx)) GO TO 100

      ind(k+1) = it
      GO TO 90
      end subroutine qsortd

!=============================================================================80
 
      subroutine qsorti(IX,ind,n)

!       Code converted using TO_F90 by Alan Miller
!       Date: 2002-12-18  Time: 11:55:47

      Integer,  INTENT(IN)    :: IX(:)
      INTEGER,  INTENT(OUT)   :: ind(:)
      INTEGER,  INTENT(IN)    :: n

!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL(wp)
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

      INTEGER     :: iu(21), il(21)
      INTEGER     :: m, i, j, k, l, ij, it, itt, indx
      REAL        :: r

      integer     :: t

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

      IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

      DO  i = 1, n
        ind(i) = i
      END DO
      m = 1
      i = 1
      j = n
      r = .375

! TOP OF LOOP

      20 IF (i >= j) GO TO 70
      IF (r <= .5898437) THEN
        r = r + .0390625
      ELSE
        r = r - .21875
      END IF

! INITIALIZE K

      30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

      ij = i + r*(j-i)
      it = ind(ij)
      t = IX(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (IX(indx) > t) THEN
        ind(ij) = indx
        ind(i) = it
        it = indx
        t = IX(it)
      END IF

! INITIALIZE L

      l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

      indx = ind(j)
      IF (IX(indx) >= t) GO TO 50
      ind(ij) = indx
      ind(j) = it
      it = indx
      t = IX(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (IX(indx) <= t) GO TO 50
      ind(ij) = indx
      ind(i) = it
      it = indx
      t = IX(it)
      GO TO 50

! INTERCHANGE ELEMENTS K AND L

      40 itt = ind(l)
      ind(l) = ind(k)
      ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

      50 l = l - 1
      indx = ind(l)
      IF (IX(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

      60 k = k + 1
      indx = ind(k)
      IF (IX(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

      IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

      IF (l-i > j-k) THEN
        il(m) = i
        iu(m) = l
        i = k
        m = m + 1
        GO TO 80
      END IF

      il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

      70 m = m - 1
      IF (m == 0) RETURN
      i = il(m)
      j = iu(m)

      80 IF (j-i >= 11) GO TO 30
      IF (i == 1) GO TO 20
      i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

      90 i = i + 1
      IF (i == j) GO TO 70
      indx = ind(i+1)
      t = IX(indx)
      it = indx
      indx = ind(i)
      IF (IX(indx) <= t) GO TO 90
      k = i

      100 ind(k+1) = ind(k)
      k = k - 1
      indx = ind(k)
      IF (t < IX(indx)) GO TO 100

      ind(k+1) = it
      GO TO 90
      end subroutine qsorti

!=============================================================================80

  subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
   
  !*****************************************************************************80
  !
  !! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
  !
  !  Discussion:
  !
  !    This is essentially a transposition operation.  
  !
  !    It is NOT an in-place algorithm.
  !
  !    This routine transposes a matrix stored in a, ja, ia format.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input, integer ( kind = 4 ) JOB, indicates whether or not to fill the values of the
  !    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
  !
  ! ipos  = starting position in ao, jao of the transposed matrix.
  !         the iao array takes this into account (thus iao(1) is set to ipos.)
  !         Note: this may be useful if one needs to append the data structure
  !         of the transpose to that of A. In this case use
  !                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
  !        for any other normal usage, enter ipos=1.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  !    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
  !    Compressed Sparse Column format.
  !
    implicit none

    integer :: n

    real (wp) a(*)
    real (wp) ao(*)
    integer :: i
    integer :: ia(n+1)
    integer :: iao(n+1)
    integer :: ipos
    integer :: j
    integer :: ja(*)
    integer :: jao(*)
    integer :: job
    integer :: k
    integer :: next
  !
  !  Compute lengths of rows of A'.
  !
    iao(1:n+1) = 0

    do i = 1, n
      do k = ia(i), ia(i+1)-1
        j = ja(k) + 1
        iao(j) = iao(j) + 1
      end do
    end do
  !
  !  Compute pointers from lengths.
  !
    iao(1) = ipos
    do i = 1, n
      iao(i+1) = iao(i) + iao(i+1)
    end do
  !
  !  Do the actual copying.
  !
    do i = 1, n
      do k = ia(i), ia(i+1)-1
        j = ja(k)
        next = iao(j)
        if ( job == 1 ) then
          ao(next) = a(k)
        end if
        jao(next) = i
        iao(j) = next + 1
      end do
    end do
  !
  !  Reshift IAO and leave.
  !
    do i = n, 1, -1
      iao(i+1) = iao(i)
    end do
    iao(1) = ipos

    return
  end subroutine csrcsc

  !============================================================================

  !============================================================================
  
  subroutine csort_block ( block_size, n, a, ja, ia, iwork, values )

  !*****************************************************************************80
  !
  !! CSORT_BLOCK sorts the elements of a CSR matrix.
  !
  !  Discussion:
  !
  !    This routine sorts the elements of a CSR matrix (stored in Compressed
  !    Sparse Row Format) in increasing order of their column indices within
  !    each row. It uses a form of bucket sort with a cost of O(nnz) where
  !    nnz = number of nonzero elements.
  !
  !    Requires an integer ( kind = 4 ) work array of size length 2*nnz.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! iwork = integer ( kind = 4 ) work array of length max ( n+1, 2*nnz )
  !         where nnz = 2* (ia(n+1)-ia(1))  ) .
  !
  ! values= logical indicating whether or not the real values a(*) must
  !         also be permuted. if (.not. values) then the array a is not
  !         touched by csort and can be a dummy array.
  !
  ! on return:
  !
  ! the matrix stored in the structure a, ja, ia is permuted in such a
  ! way that the column indices are in increasing order within each row.
  ! iwork(1:nnz) contains the permutation used  to rearrange the elements.
  !
    implicit none

    integer ::  n

    integer :: block_size
    real(wp) :: a(block_size,block_size,*)
    integer :: i
    integer :: ia(n+1)
    integer :: ifirst
    integer :: irow
    integer :: iwork(*)
    integer :: j
    integer :: ja(*)
    integer :: k
    integer :: ko
    integer :: next
    integer :: nnz
    logical :: values
  !
  !  Count the number of elements in each column.
  !
    iwork(1:n+1) = 0

    do i = 1, n
       do k = ia(i), ia(i+1)-1
          j = ja(k) + 1
          iwork(j) = iwork(j) + 1
       end do
    end do
  !
  !  Compute pointers from lengths.
  !
    iwork(1) = 1

    do i = 1, n
       iwork(i+1) = iwork(i) + iwork(i+1)
    end do
  !
  !  Get the positions of the nonzero elements in order of columns.
  !
    ifirst = ia(1)
    nnz = ia(n+1)-ifirst

    do i = 1, n
      do k = ia(i), ia(i+1)-1
        j = ja(k)
        next = iwork(j)
        iwork(nnz+next) = k
        iwork(j) = next + 1
      end do
    end do
  !
  !  Convert to coordinate format.
  !
    do i = 1, n
      do k = ia(i), ia(i+1)-1
        iwork(k) = i
      end do
    end do
  !
  !  Loop to find permutation: for each element find the correct
  !  position in (sorted) arrays A, JA.  Record this in IWORK.
  !
    do k = 1, nnz
       ko = iwork(nnz+k)
       irow = iwork(ko)
       next = ia(irow)
  !
  !  The current element should go in next position in row. IWORK
  !  records this position.
  !
       iwork(ko) = next
       ia(irow) = next + 1
    end do
  !
  !  Perform an in-place permutation of the arrays.
  !
       call ivperm ( nnz, ja(ifirst), iwork )

       if ( values ) then
         call dvperm_block ( block_size, nnz, a(:,:,ifirst), iwork )
       end if
  !
  !  Reshift the pointers of the original matrix back.
  !
    do i = n, 1, -1
      ia(i+1) = ia(i)
    end do

    ia(1) = ifirst

    return
  end subroutine csort_block

  !============================================================================

  !============================================================================
  
  subroutine dvperm_block ( block_size, n, x, perm )

  !*****************************************************************************80
  !
  !! DVPERM performs an in-place permutation of a real vector.
  !
  !  Discussion:
  !
  !    This routine permutes a real vector X using a permutation PERM.
  !
  !    On return, the vector X satisfies,
  !
  !      x(perm(j)) :== x(j), j = 1,2,.., n
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the length of X.
  !
  !    Input/output, real X(N), the vector to be permuted.
  !
  !    Input, integer ( kind = 4 ) PERM(N), the permutation.
  !
    implicit none

    integer :: block_size
    integer :: n

    integer :: ii
    integer :: init
    integer :: k
    integer :: next
    integer :: perm(n)
    real(wp), dimension(block_size,block_size) :: tmp
    real(wp), dimension(block_size,block_size) :: tmp1
    real(wp), dimension(block_size,block_size,n) :: x

    init = 1
    tmp = x(:,:,init)
    ii = perm(init)
    perm(init)= -perm(init)
    k = 0
  !
  !  The main loop.
  !
   6  continue

     k = k + 1
  !
  !  Save the chased element.
  !
    tmp1 = x(:,:,ii)
    x(:,:,ii) = tmp
    next = perm(ii)

    if ( next < 0 ) then
      go to 65
    end if
  !
  !  Test for end.
  !
    if ( n < k ) then
      perm(1:n) = -perm(1:n)
      return
    end if

    tmp = tmp1
    perm(ii) = -perm(ii)
    ii = next
  !
  !  End of the loop.
  !
    go to 6
  !
  !  Reinitialize cycle.
  !
   65   continue

    init = init + 1

    if ( n < init ) then 
      perm(1:n) = -perm(1:n)
      return
    end if

    if ( perm(init) < 0 ) then
      go to 65
    end if

    tmp = x(:,:,init)
    ii = perm(init)
    perm(init) = -perm(init)
    go to 6

  end subroutine dvperm_block

  !============================================================================

  subroutine a_bl_mu_x(block_size,n_row_blocks,x,y,a,ja,ia) 
    integer :: n_row_blocks, block_size
    real(wp) :: x(block_size*n_row_blocks),y(block_size*n_row_blocks), &
      & a(block_size,block_size,*) 
    integer :: ja(*), ia(*)

    !-----------------------------------------------------------------------
    !         A block sparse row times a vector
    !----------------------------------------------------------------------- 
    ! multiplies a matrix by a vector using the dot product form
    ! Matrix A is stored in block compressed sparse row storage.
    !
    ! on entry:
    !----------
    ! n_row_blocks = row block dimension of A
    ! x     = real array of length equal to the block column dimension of
    !         the A matrix.
    ! a, ja,
    !    ia = input matrix in block compressed sparse row format.
    !
    ! on return:
    !-----------
    ! y     = real array of length n, containing the product y=Ax
    !
    !-----------------------------------------------------------------------
    
    real(wp) :: tmp
    integer :: i_row_block, i_column_block, i_row, i_column, cnt 

    cnt = 1

    do i_row_block = 1, n_row_blocks
      do i_column_block = ia(i_row_block), ia(i_row_block+1)-1
        do i_row = 1, block_size
          tmp = 0.0
          do i_column = 1, block_size
            tmp = tmp + a(i_row,i_column,cnt)* &
              & x((ja(i_column_block)-1)*block_size + i_column)
          enddo
          y(i_row + (i_row_block-1)*block_size) = &
            & y(i_row + (i_row_block-1)*block_size) + tmp
        enddo
        cnt = cnt + 1
      enddo
    enddo
    
 
  end subroutine a_bl_mu_x
  
  !============================================================================


end module unary_mod
