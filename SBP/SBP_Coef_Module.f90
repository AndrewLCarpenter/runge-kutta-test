      module SBP_Coef_Module

      implicit none

      private
      public  ::  D1_242

      contains

      subroutine D1_242(n,nnz,ia,ja,a,Pmat,Pinv)

      implicit none

      integer,   parameter                      :: wp = 8

      integer,                    intent(in   ) :: n,nnz

      integer,   dimension(n+1),  intent(  out) :: ia
      integer,   dimension(nnz),  intent(  out) :: ja

      real(wp),  dimension(nnz),  intent(  out) ::  a
      real(wp),  dimension(n  ),  intent(  out) ::  Pmat, Pinv

      integer                                   :: i,j
      integer                                   :: icnt, jcnt
      real(wp)                                  :: h
      real(wp),  allocatable, dimension(:)      :: dmat
      real(wp),  allocatable, dimension(:,:)    :: D1blk, D2blk
      real(wp),  parameter                      :: tol=1.0e-12_wp

      continue

      allocate(dmat(1:5))
      allocate(D1blk(1:4,1:6))
      allocate(D2blk(1:4,1:6))

      h = 1.0_wp / (n - 1)

      Pmat(1:  4) = reshape(                              &
                  & (/17.0_wp/48.0_wp, 59.0_wp/48.0_wp,   &
                  &   43.0_wp/48.0_wp, 49.0_wp/48.0_wp /),&
                  & (/4/)  )
      Pmat(5:n-4) = 1.0_wp
      Pmat(n-3:n) = reshape(                              &
                    (/49.0_wp/48.0_wp, 43.0_wp/48.0_wp,   &
                      59.0_wp/48.0_wp, 17.0_wp/48.0_wp /),&
                  & (/4/)  )

      Pinv(:) = 1.0_wp / Pmat(:)

      dmat(:) = reshape(                                               &
              & (/ 1.0_wp/12,-8.0_wp/12,0.0_wp,8.0_wp/12,-1.0_wp/12 /),&
              & (/5/))

      D1blk   = reshape(                                                                 &
              & (/-24.0_wp/17, 59.0_wp/34, -4.0_wp/17,-3.0_wp/34, 0.0_wp   , 0.0_wp,     &
              &    -1.0_wp/2 ,  0.0_wp   ,  1.0_wp/2 , 0.0_wp   , 0.0_wp   , 0.0_wp,     &
              &     4.0_wp/43,-59.0_wp/86,  0.0_wp   ,59.0_wp/86,-4.0_wp/43, 0.0_wp,     &
              &     3.0_wp/98,  0.0_wp   ,-59.0_wp/98, 0.0_wp   ,32.0_wp/49,-4.0_wp/49/),&
              & (/4,6/)  )

      do i=1,4
         do j=1,6
           D2blk(i,j)= -D1blk(5-i,7-j)
         end do
      end do

      !  initialize CSR matrix information
      icnt  = 0      ; jcnt = 0 

      ia(:) = 0      ; ja(:) = 0      ;
       a(:) = 0.0_wp ;
      ia(1) = 1      ! start at beginning of array

      !  load derivative matrix info into CSR buckets:  Left block boundary

      do i = 1,4
        jcnt = 0 
        do j = 1,6
          if(abs(D1blk(i,j)) >= tol) then
              icnt = icnt + 1 ; jcnt = jcnt + 1 ;
          ja(icnt) = j
           a(icnt) = D1blk(i,j)
          endif
        enddo
        ia(i+1) = ia(i) + jcnt
      enddo
  
      !  load derivative matrix info into CSR buckets:  interior
      do i = 5,n-4
        jcnt = 0 
        do j = 1,5
          if(abs(dmat(j)) >= tol) then
              icnt = icnt + 1 ; jcnt = jcnt + 1 ;
          ja(icnt) = n-6+j
           a(icnt) = dmat(j)
          endif
        enddo
        ia(i+1) = ia(i) + jcnt
      enddo
 
      do i = 1,4
        jcnt = 0 
        do j = 1,6
          if(abs(D2blk(i,j)) >= tol) then
              icnt = icnt + 1 ; jcnt = jcnt + 1 ;
          ja(icnt) = n-6+j
           a(icnt) = D2blk(i,j)
          endif
        enddo
        ia(i+1) = ia(i) + jcnt
      enddo

      a(:) = a(:) / h

      deallocate(dmat,D1blk,D2blk)

      return

      end subroutine D1_242

      end module SBP_Coef_Module
