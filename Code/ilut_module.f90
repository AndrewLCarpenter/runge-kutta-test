      module ilut_module

      use precision_vars, only : wp
      use matvec_module , only : amux

      implicit none; save

      private

      public ::   ilut, ilutp, ilu0, iluk, lusol, ilud, iludp
      public ::   ilu0_csr, lusol_csr, msrcsr
      public ::   lusol_M_rhs

      contains

!=============================================================================80

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                   ITERATIVE SOLVERS MODULE                           c
!----------------------------------------------------------------------c
! This Version Dated: August 13, 1996. Warning: meaning of some        c
! ============ arguments have changed w.r.t. earlier versions. Some    c
!              Calling sequences may also have changed                 c
!----------------------------------------------------------------------c
! Contents:                                                            c
!-------------------------preconditioners------------------------------c
!                                                                      c
! ILUT    : Incomplete LU factorization with dual truncation strategy  c
! ILUTP   : ILUT with column  pivoting                                 c
! ILUD    : ILU with single dropping + diagonal compensation (~MILUT)  c
! ILUDP   : ILUD with column pivoting                                  c
! ILUK    : level-k ILU                                                c
! ILU0    : simple ILU(0) preconditioning                              c
! MILU0   : MILU(0) preconditioning                                    c
!                                                                      c
!----------sample-accelerator-and-LU-solvers---------------------------c
!                                                                      c
! PGMRES  : preconditioned GMRES solver                                c
! LUSOL   : forward followed by backward triangular solve (Precond.)   c
! LUTSOL  : solving v = (LU)^{-T} u (used for preconditioning)         c
!                                                                      c
!-------------------------utility-routine------------------------------c
!                                                                      c
! QSPLIT  : quick split routine used by ilut to sort out the k largest c
!           elements in absolute value                                 c
!                                                                      c
!----------------------------------------------------------------------c
!                                                                      c
! Note: all preconditioners are preprocessors to pgmres.               c
! usage: call preconditioner then call pgmres                          c
!                                                                      c
!----------------------------------------------------------------------c
      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr) 
!-----------------------------------------------------------------------
      implicit none 

      integer,                  intent(in   ) :: n, lfil, iwk
      real(wp),                 intent(in   ) :: droptol 
      real(wp), dimension(*),   intent(in   ) :: a
      real(wp), dimension(*),   intent(inout) :: alu
      real(wp), dimension(n+1), intent(inout) :: w

      integer,     dimension(*),   intent(in   ) :: ja
      integer,     dimension(n+1), intent(in   ) :: ia
      integer,     dimension(*),   intent(inout) :: jlu
      integer,     dimension(n),   intent(inout) :: ju
      integer,     dimension(2*n), intent(inout) :: jw
      integer,                     intent(  out) :: ierr 

!----------------------------------------------------------------------*
!                      *** ILUT preconditioner ***                     *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
!     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!----------------------------------------------------------------------*
! PARAMETERS                                                            
!-----------                                                            
!                                                                       
! on entry:                                                             
!==========                                                             
! n       = integer. The row dimension of the matrix A. The matrix      
!                                                                       
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!                                                                       
! lfil    = integer. The fill-in parameter. Each row of L and each row  
!           of U will have a maximum of lfil elements (excluding the    
!           diagonal element). lfil must be  >=  0.                     
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO 
!           EARLIER VERSIONS.                                           
!                                                                       
! droptol = real(wp). Sets the threshold for dropping small terms in the  
!           factorization. See below for details on dropping strategy.  
!                                                                       
!                                                                       
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays   
!           are not big enough to store the ILU factorizations, ilut    
!           will stop with an error message.                            
!                                                                       
! On return:                                                            
!===========                                                            
!                                                                       
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju      = integer array of length n containing the pointers to        
!           the beginning of each row of U in the matrix alu,jlu.       
!                                                                       
! ierr    = integer. Error message with the following meaning.          
!           ierr  =  0   --> successful return.                         
!           ierr  >  0   --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.          
!                            (The elimination process has generated a   
!                            row in L or U whose length is  >   n.)    
!           ierr  = -2   --> The matrix L overflows the array al.       
!           ierr  = -3   --> The matrix U overflows the array alu.      
!           ierr  = -4   --> Illegal value for lfil.                    
!           ierr  = -5   --> zero row encountered.                      
!                                                                       
! work arrays:                                                          
!=============                                                          
! jw      = integer work array of length 2*n.                           
! w       = real work array of length n+1.                              
!                                                                       
!---------------------------------------------------------------------- 
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]       
! jw(n+1:2n)  stores nonzero indicators                                 
!                                                                       
! Notes:                                                                
! ------                                                                
! The diagonal elements of the input matrix must be  nonzero (at least  
! 'structurally').                                                      
!                                                                       
!----------------------------------------------------------------------*
!---- Dual drop strategy works as follows.                             *
!                                                                      *
!     1) Theresholding in L and U as set by droptol. Any element whose *
!        magnitude is less than some tolerance (relative to the abs    *
!        value of diagonal element in u) is dropped.                   *
!                                                                      *
!     2) Keeping only the largest lfil elements in the i-th row of L   *
!        and the largest lfil elements in the i-th row of U (excluding *
!        diagonal elements).                                           *
!                                                                      *
! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
! keeping  the largest  elements in  each row  of L  and U.   Taking   *
! droptol  /=   0 but lfil=n will give  the usual threshold strategy   *
! (however, fill-in is then mpredictible).                             *
!----------------------------------------------------------------------*
!     locals                                                            
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,lngth 
      integer icomp
      real(wp) tnorm, t, abs, s, fact , eps

      icomp = 0

      if (lfil  <  0) goto 998 
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)    
!     and pointer array.                                                
!-----------------------------------------------------------------------
      ju0    = n+2 
      jlu(1) = ju0 
                                                                        
!     initialize nonzero indicator array.                               
                                                                        
      do j = 1,n 
         jw(n+j)  = 0 
      enddo

!-----------------------------------------------------------------------
      do ii = 1, n          !     beginning of main loop.                                           
!-----------------------------------------------------------------------
         j1 = ia(ii  ) 
         j2 = ia(ii+1) - 1 
         tnorm = 0.0_wp 
         do k=j1,j2 
            tnorm = tnorm+abs(a(k)) 
         enddo
         if (tnorm  ==  0.0_wp) goto 999 
         tnorm = tnorm/real(j2-j1+1) 
                                                                        
!        unpack L-part and U-part of row of A in arrays w                  
                                                                        
         lenu = 1 ; lenl = 0 ;
         jw(  ii) = ii ;
         jw(n+ii) = ii ;

         w(ii) = 0.0 

         do j = j1, j2 
            k = ja(j) 
            if(k == ii) then
              t = a(j) + 1.0e-05_wp
            else
              t = a(j)
            endif
            if (k  <  ii) then 
               lenl = lenl+1 
               jw(lenl) = k 
                w(lenl) = t 
               jw(n+k)  = lenl 
            else if (k  ==  ii) then 
               w(ii) = t 
            else 
               lenu = lenu+1 
               jpos = ii+lenu-1 
               jw(jpos) = k 
                w(jpos) = t 
               jw(n+k)  = jpos 
            endif 
         enddo

         jj    = 0 ; lngth = 0  ;

!        eliminate previous rows                                           

  150    jj = jj+1 
         if (jj  >  lenl) goto 160 
!-----------------------------------------------------------------------
!        To do the elimination in the correct order we must select
!        the smallest column index among jw(k), k=jj+1, ..., lenl.         
!-----------------------------------------------------------------------
         jrow = jw(jj) 
            k = jj 
                                                                        
!        determine smallest column index                                   
                                                                        
         do j=jj+1,lenl 
            if (jw(j)  <  jrow) then 
               jrow = jw(j) 
               k = j 
            endif 
         enddo 
                                                                        
         if (k  /=  jj) then 
!           exchange in jw                                                    
            j      = jw(jj) 
            jw(jj) = jw(k) 
            jw(k)  = j 
!           exchange in jr                                                    
            jw(n+jrow) = jj 
            jw(n+j)    = k 
!           exchange in w                                                     
            s = w(jj) 
            w(jj) = w(k) 
            w(k) = s 
         endif 
                                                                        
!        zero out element in row by setting jw(n+jrow) to zero.            
                                                                        
         jw(n+jrow) = 0 
                                                                        
!        get the multiplier for row to be eliminated (jrow).               
                                                                        
         fact = w(jj)*alu(jrow) 
         if (abs(fact)  <=  droptol) goto 150 
                                                                        
!        combine current row and row jrow                                  
                                                                        
         do k = ju(jrow), jlu(jrow+1)-1 
            s = fact*alu(k) 
            j = jlu(k) 
            jpos = jw(n+j) 
            if (j  >=  ii) then    !  dealing with upper part.                                          
               if (jpos  ==  0) then  ! this is a fill-in element 
                  lenu = lenu+1 
                  if (lenu  >  n) goto 995 
                  i = ii+lenu-1 
                  jw(i)   = j 
                  jw(n+j) = i 
                  w(i) = - s 
               else                   ! this is not a fill-in element
                  w(jpos) = w(jpos) - s
               endif 
            else                  !  dealing  with lower part.                                         
               if (jpos  ==  0) then  ! this is a fill-in element
                  lenl = lenl+1 
                  if (lenl  >  n) goto 995 
                  jw(lenl) = j 
                  jw(n+j) = lenl 
                  w(lenl) = - s 
               else                   ! this is not a fill-in element
                  w(jpos) = w(jpos) - s 
               endif 
            endif 
         enddo

!        store this pivot element -- (from left to right -- no danger of   
!        overlap with the working elements in L (pivots).                  
                                                                        
         lngth = lngth+1 
         w(lngth) = fact 
         jw(lngth)  = jrow 
         goto 150 
  160    continue 
                                                                        
!        reset double-pointer to zero (U-part)                             
                                                                        
         do k = 1, lenu 
            jw(n+jw(ii+k-1)) = 0 
         enddo

!        update L-matrix                                                   
                                                                        
         lenl  = lngth
         lngth = min0(lenl,lfil) 
                                                                        
!        sort by quick-split                                               
                                                                        
         call qsplit (w,jw,lenl,lngth) 
                                                                        
!        store L-part                                                      
                                                                        
         do k = 1, lngth
            if (ju0  >  iwk) goto 996 
            alu(ju0) =  w(k) 
            jlu(ju0) =  jw(k) 
            ju0 = ju0+1 
         enddo
                                                                        
!        save pointer to beginning of row ii of U                          
                                                                        
         ju(ii) = ju0 
                                                                        
!        update U-matrix -- first apply dropping strategy                  
                                                                        
         lngth = 0 
         do k = 1, lenu-1 
            if (abs(w(ii+k))  >  droptol*tnorm) then 
               lngth = lngth+1 
               w(ii+lngth) = w(ii+k) 
               jw(ii+lngth) = jw(ii+k) 
            endif 
         enddo 
         lenu = lngth+1 
         lngth = min0(lenu,lfil) 
 
         call qsplit (w(ii+1), jw(ii+1), lenu-1,lngth) 
                                                                        
!-----------------------------------------------------------------------

         if (lngth + ju0 >  iwk) goto 997

!        store U-part
                                                                        
!        t = abs(w(ii)) 
         do k = ii+1,ii+lngth-1 
            jlu(ju0) = jw(k) 
            alu(ju0) = w(k) 
!           t = t + abs(w(k)) 
            ju0 = ju0+1 
         enddo
                                                                        
!        store inverse of diagonal element of u                            

         eps = 1.0e-13_wp
         if(abs(w(ii)) < eps*tnorm)then
           icomp = icomp + 1
           if(w(ii) > 0.0_wp ) w(ii) = w(ii) + eps*tnorm
           if(w(ii) < 0.0_wp ) w(ii) = w(ii) - eps*tnorm
         endif
                                                                        
         alu(ii) = 1.0_wp/ w(ii) 
                                                                        
!        update pointer to beginning of next row of U.                     
                                                                        
         jlu(ii+1) = ju0 

!-----------------------------------------------------------------------
      enddo                 !     ending of main loop
!-----------------------------------------------------------------------

      if(icomp /= 0)write(*,*)'diagonal compensation', icomp

      ierr =  0  !  Normal Return
      return 
  995 ierr = -1    ! incomprehensible error. Matrix must be wrong.                     
      return 
  996 ierr = -2    ! insufficient storage in L.                                        
      return 
  997 ierr = -3    ! insufficient storage in U.                                        
      return 
  998 ierr = -4    ! illegal lfil entered.                                             
      return 
  999 ierr = -5    ! zero row encountered                                              
      return 

      end subroutine ilut                                           

!=============================================================================80

!---------------------------------------------------------------------- 
      subroutine ilutp(n,a,ja,ia,lfil,droptol,permtol,mbloc,alu,        &
     &     jlu,ju,iwk,w,jw,iperm,ierr)                                  
!-----------------------------------------------------------------------
!     implicit none                                                     
      integer n,ja(*),ia(n+1),lfil,jlu(*),ju(n),jw(2*n),iwk,            &
     &     iperm(2*n),ierr                                              
      real(wp) a(*), alu(*), w(n+1), droptol 

!----------------------------------------------------------------------*
!       *** ILUTP preconditioner -- ILUT with pivoting  ***            *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
! author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996.     *
!----------------------------------------------------------------------*
! on entry:                                                             
!==========                                                             
! n       = integer. The dimension of the matrix A.                     
!                                                                       
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!           ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR      
!           DETAILS.                                                    
!                                                                       
! lfil    = integer. The fill-in parameter. Each row of L and each row  
!           of U will have a maximum of lfil elements (excluding the    
!           diagonal element). lfil must be  >=  0.                     
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO 
!           EARLIER VERSIONS.                                           
!                                                                       
! droptol = real(wp). Sets the threshold for dropping small terms in the  
!           factorization. See below for details on dropping strategy.  
!                                                                       
! permtol = tolerance ratio used to  determne whether or not to permute 
!           two columns.  At step i columns i and j are permuted when   
!                                                                       
!                     abs(a(i,j))*permtol  >  abs(a(i,i))              
!                                                                       
!           [0 --> never permute; good values 0.1 to 0.01]              
!                                                                       
! mbloc   = if desired, permuting can be done only within the diagonal  
!           blocks of size mbloc. Useful for PDE problems with several  
!           degrees of freedom.. If feature not wanted take mbloc=n.    
!                                                                       
!                                                                       
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays   
!           are not big enough to store the ILU factorizations, ilut    
!           will stop with an error message.                            
!                                                                       
! On return:                                                            
!===========                                                            
!                                                                       
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju      = integer array of length n containing the pointers to        
!           the beginning of each row of U in the matrix alu,jlu.       
!                                                                       
! iperm   = contains the permutation arrays.                            
!           iperm(1:n) = old numbers of unknowns                        
!           iperm(n+1:2*n) = reverse permutation = new unknowns.        
!                                                                       
! ierr    = integer. Error message with the following meaning.          
!           ierr  =  0   --> successful return.                         
!           ierr  >  0   --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.          
!                            (The elimination process has generated a   
!                            row in L or U whose length is  >   n.)    
!           ierr  = -2   --> The matrix L overflows the array al.       
!           ierr  = -3   --> The matrix U overflows the array alu.      
!           ierr  = -4   --> Illegal value for lfil.                    
!           ierr  = -5   --> zero row encountered.                      
!                                                                       
! work arrays:                                                          
!=============                                                          
! jw      = integer work array of length 2*n.                           
! w       = real work array of length n                                 
!                                                                       
! IMPORTANR NOTE:                                                       
! --------------                                                        
! TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE,     
! THE MATRIX A IS PERMUTED ON RETURN. [all column indices are           
! changed]. SIMILARLY FOR THE U MATRIX.                                 
! To permute the matrix back to its original state use the loop:        
!                                                                       
!      do k=ia(1), ia(n+1)-1                                            
!         ja(k) = iperm(ja(k))                                          
!      enddo                                                            
!                                                                       
!-----------------------------------------------------------------------
!     local variables                                                   
!                                                                       
      integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,mbloc, &
     &     icut                                                         
      real(wp) s, tmp, tnorm,xmax,xmax0, fact, abs, t, permtol 
!                                                                       
      if (lfil  <  0) goto 998 
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)    
!     and pointer array.                                                
!-----------------------------------------------------------------------
      ju0 = n+2 
      jlu(1) = ju0 
!                                                                       
!  integer double pointer array.                                        
!                                                                       
      do 1 j=1, n 
         jw(n+j)    = 0 
         iperm(  j) = j 
         iperm(n+j) = j 
    1 continue 
!-----------------------------------------------------------------------
!     beginning of main loop.                                           
!-----------------------------------------------------------------------
      do 500 ii = 1, n 
         j1 = ia(ii) 
         j2 = ia(ii+1) - 1 
         tnorm = 0.0_wp 
         do 501 k=j1,j2 
            tnorm = tnorm+abs(a(k)) 
  501    continue 
         if (tnorm  ==  0.0_wp) goto 999 
         tnorm = tnorm/(j2-j1+1) 
!                                                                       
!     unpack L-part and U-part of row of A in arrays  w  --             
!                                                                       
         lenu = 1 
         lenl = 0 
         jw(  ii) = ii 
          w(  ii) = 0.0_wp 
         jw(n+ii) = ii 
!                                                                       
         do 170  j = j1, j2 
            k = iperm(n+ja(j)) 
            t = a(j) 
            if (k  <  ii) then 
               lenl = lenl+1 
               jw(lenl) = k 
                w(lenl) = t 
               jw(n+k)  = lenl 
            else if (k  ==  ii) then 
               w(ii) = t 
            else 
               lenu = lenu+1 
               jpos = ii+lenu-1 
               jw(jpos) = k 
                w(jpos) = t 
               jw(n+k)  = jpos 
            endif 
  170    continue 
         jj = 0 
         len = 0 
!                                                                       
!     eliminate previous rows                                           
!                                                                       
  150    jj = jj+1 
         if (jj  >  lenl) goto 160 
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.         
!-----------------------------------------------------------------------
         jrow = jw(jj) 
         k = jj 
!                                                                       
!     determine smallest column index                                   
!                                                                       
         do 151 j=jj+1,lenl 
            if (jw(j)  <  jrow) then 
               jrow = jw(j) 
               k = j 
            endif 
  151    continue 
!                                                                       
         if (k  /=  jj) then 
!     exchange in jw                                                    
            j = jw(jj) 
            jw(jj) = jw(k) 
            jw(k) = j 
!     exchange in jr                                                    
            jw(n+jrow) = jj 
            jw(n+j) = k 
!     exchange in w                                                     
            s = w(jj) 
            w(jj) = w(k) 
            w(k) = s 
         endif 
!                                                                       
!     zero out element in row by resetting jw(n+jrow) to zero.          
!                                                                       
         jw(n+jrow) = 0 
!                                                                       
!     get the multiplier for row to be eliminated: jrow                 
!                                                                       
         fact = w(jj)*alu(jrow) 
!                                                                       
!     drop term if small                                                
!                                                                       
         if (abs(fact)  <=  droptol) goto 150 
!                                                                       
!     combine current row and row jrow                                  
!                                                                       
         do 203 k = ju(jrow), jlu(jrow+1)-1 
            s = fact*alu(k) 
!     new column number                                                 
            j = iperm(n+jlu(k)) 
            jpos = jw(n+j) 
            if (j  >=  ii) then 
!                                                                       
!     dealing with upper part.                                          
!                                                                       
               if (jpos  ==  0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenu = lenu+1 
                  i = ii+lenu-1 
                  if (lenu  >  n) goto 995 
                  jw(i) = j 
                  jw(n+j) = i 
                  w(i) = - s 
               else 
!     no fill-in element --                                             
                  w(jpos) = w(jpos) - s 
               endif 
            else 
!                                                                       
!     dealing with lower part.                                          
!                                                                       
               if (jpos  ==  0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                 lenl = lenl+1 
                 if (lenl  >  n) goto 995 
                 jw(lenl) = j 
                 jw(n+j) = lenl 
                 w(lenl) = - s 
              else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                 w(jpos) = w(jpos) - s 
              endif 
           endif 
  203       continue 
!                                                                       
!     store this pivot element -- (from left to right -- no danger of   
!     overlap with the working elements in L (pivots).                  
!                                                                       
        len = len+1 
        w(len) = fact 
        jw(len)  = jrow 
        goto 150 
  160   continue 
!                                                                       
!     reset double-pointer to zero (U-part)                             
!                                                                       
        do 308 k=1, lenu 
           jw(n+jw(ii+k-1)) = 0 
  308       continue 
!                                                                       
!     update L-matrix                                                   
!                                                                       
        lenl = len 
        len = min0(lenl,lfil) 
!                                                                       
!     sort by quick-split                                               
!                                                                       
        call qsplit (w,jw,lenl,len) 
!                                                                       
!     store L-part -- in original coordinates ..                        
!                                                                       
        do 204 k=1, len 
           if (ju0  >  iwk) goto 996 
           alu(ju0) =  w(k) 
           jlu(ju0) = iperm(jw(k)) 
           ju0 = ju0+1 
  204   continue 
!                                                                       
!     save pointer to beginning of row ii of U                          
!                                                                       
        ju(ii) = ju0 
!                                                                       
!     update U-matrix -- first apply dropping strategy                  
!                                                                       
         len = 0 
         do k=1, lenu-1 
            if (abs(w(ii+k))  >  droptol*tnorm) then 
               len = len+1 
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif 
         enddo 
         lenu = len+1 
         len = min0(lenu,lfil) 
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len) 
!                                                                       
!     determine next pivot --                                           
!                                                                       
        imax = ii 
        xmax = abs(w(imax)) 
        xmax0 = xmax 
        icut = ii - 1 + mbloc - mod(ii-1,mbloc) 
        do k=ii+1,ii+len-1 
           t = abs(w(k)) 
           if (t  >  xmax .and. t*permtol  >  xmax0 .and.             &
     &          jw(k)  <=  icut) then                                   
              imax = k 
              xmax = t 
           endif 
        enddo 
!                                                                       
!     exchange w's                                                      
!                                                                       
        tmp = w(ii) 
        w(ii) = w(imax) 
        w(imax) = tmp 
!                                                                       
!     update iperm and reverse iperm                                    
!                                                                       
        j = jw(imax) 
        i = iperm(ii) 
        iperm(ii) = iperm(j) 
        iperm(j) = i 
!                                                                       
!     reverse iperm                                                     
!                                                                       
        iperm(n+iperm(ii)) = ii 
        iperm(n+iperm(j)) = j 
!-----------------------------------------------------------------------
!                                                                       
        if (len + ju0  >  iwk) goto 997 
!                                                                       
!     copy U-part in original coordinates                               
!                                                                       
        do 302 k=ii+1,ii+len-1 
           jlu(ju0) = iperm(jw(k)) 
           alu(ju0) = w(k) 
           ju0 = ju0+1 
  302       continue 
!                                                                       
!     store inverse of diagonal element of u                            
!                                                                       
        if (w(ii)  ==  0.0) w(ii) = (1.0D-4 + droptol)*tnorm 
        alu(ii) = 1.0_wp/ w(ii) 
!                                                                       
!     update pointer to beginning of next row of U.                     
!                                                                       
        jlu(ii+1) = ju0 
!-----------------------------------------------------------------------
!     end main loop                                                     
!-----------------------------------------------------------------------
  500 continue 
                                                                        
!     permute all column indices of LU ...                              
                                                                        
      do k = jlu(1),jlu(n+1)-1 
         jlu(k) = iperm(n+jlu(k)) 
      enddo 
                                                                        
!     ...and of A                                                       
                                                                        
      do k=ia(1), ia(n+1)-1 
         ja(k) = iperm(n+ja(k)) 
      enddo 
!                                                                       
      ierr = 0 
      return 
!                                                                       
!     incomprehensible error. Matrix must be wrong.                     
!                                                                       
  995 ierr = -1 
      return 
!                                                                       
!     insufficient storage in L.                                        
!                                                                       
  996 ierr = -2 
      return 
!                                                                       
!     insufficient storage in U.                                        
!                                                                       
  997 ierr = -3 
      return 
!                                                                       
!     illegal lfil entered.                                             
!                                                                       
  998 ierr = -4 
      return 
!                                                                       
!     zero row encountered                                              
!                                                                       
  999 ierr = -5 
      return 

      end subroutine ilutp                                           

!=============================================================================80

!-----------------------------------------------------------------------
      subroutine ilud(n,a,ja,ia,alph,tol,alu,jlu,ju,iwk,w,jw,ierr) 
!-----------------------------------------------------------------------
      implicit none 
      integer n 
      real(wp) a(*),alu(*),w(2*n),tol, alph 
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),iwk,ierr 
!----------------------------------------------------------------------*
!                     *** ILUD preconditioner ***                      *
!    incomplete LU factorization with standard droppoing strategy      *
!----------------------------------------------------------------------*
! Author: Yousef Saad * Aug. 1995 --                                   *
!----------------------------------------------------------------------*
! This routine computes the ILU factorization with standard threshold  *
! dropping: at i-th step of elimination, an element a(i,j) in row i is *
! dropped  if it satisfies the criterion:                              *
!                                                                      *
!  abs(a(i,j)) < tol * [average magnitude of elements in row i of A]   *
!                                                                      *
! There is no control on memory size required for the factors as is    *
! done in ILUT. This routines computes also various diagonal compensa- *
! tion ILU's such MILU. These are defined through the parameter alph   *
!----------------------------------------------------------------------*
! on entry:                                                             
!==========                                                             
! n       = integer. The row dimension of the matrix A. The matrix      
!                                                                       
! a,ja,ia = matrix stored in Compressed Sparse Row format               
!                                                                       
! alph    = diagonal compensation parameter -- the term:                
!                                                                       
!           alph*(sum of all dropped out elements in a given row)       
!                                                                       
!           is added to the diagonal element of U of the factorization  
!           Thus: alph = 0 ---> ~ ILU with threshold,                   
!                 alph = 1 ---> ~ MILU with threshold.                  
!                                                                       
! tol     = Threshold parameter for dropping small terms in the         
!           factorization. During the elimination, a term a(i,j) is     
!           dropped whenever abs(a(i,j))  <  tol * [weighted norm of   
!           row i]. Here weighted norm = 1-norm / number of nnz         
!           elements in the row.                                        
!                                                                       
! iwk     = The length of arrays alu and jlu -- this routine will stop  
!           if storage for the factors L and U is not sufficient        
!                                                                       
! On return:                                                            
!===========                                                            
!                                                                       
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju      = integer array of length n containing the pointers to        
!           the beginning of each row of U in the matrix alu,jlu.       
!                                                                       
! ierr    = integer. Error message with the following meaning.          
!           ierr  = 0    --> successful return.                         
!           ierr  >  0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.          
!                            (The elimination process has generated a   
!                            row in L or U whose length is  >   n.)    
!           ierr  = -2   --> Insufficient storage for the LU factors -- 
!                            arrays alu/ jalu are  overflowed.          
!           ierr  = -3   --> Zero row encountered.                      
!                                                                       
! Work Arrays:                                                          
!=============                                                          
! jw      = integer work array of length 2*n.                           
! w       = real work array of length n                                 
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]       
! jw(n+1:2n)  stores the nonzero indicator.                             
!                                                                       
! Notes:                                                                
! ------                                                                
! All diagonal elements of the input matrix must be  nonzero.           
!                                                                       
!-----------------------------------------------------------------------
!     locals                                                            
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      real(wp) tnorm, t, abs, s, fact, dropsum 
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)    
!     and pointer array.                                                
!-----------------------------------------------------------------------
      ju0 = n+2 
      jlu(1) = ju0 
!                                                                       
!     initialize nonzero indicator array.                               
!                                                                       
      do 1 j=1,n 
         jw(n+j)  = 0 
    1 continue 
!-----------------------------------------------------------------------
!     beginning of main loop.                                           
!-----------------------------------------------------------------------
      do 500 ii = 1, n 
         j1 = ia(ii) 
         j2 = ia(ii+1) - 1 
         dropsum = 0.0_wp 
         tnorm = 0.0_wp 
         do 501 k=j1,j2 
            tnorm = tnorm + abs(a(k)) 
  501    continue 
         if (tnorm  ==  0.0) goto 997 
         tnorm = tnorm / real(j2-j1+1) 
!                                                                       
!     unpack L-part and U-part of row of A in arrays w                  
!                                                                       
         lenu = 1 
         lenl = 0 
         jw(ii) = ii 
         w(ii) = 0.0 
         jw(n+ii) = ii 
!                                                                       
         do 170  j = j1, j2 
            k = ja(j) 
            t = a(j) 
            if (k  <  ii) then 
               lenl = lenl+1 
               jw(lenl) = k 
               w(lenl) = t 
               jw(n+k) = lenl 
            else if (k  ==  ii) then 
               w(ii) = t 
            else 
               lenu = lenu+1 
               jpos = ii+lenu-1 
               jw(jpos) = k 
               w(jpos) = t 
               jw(n+k) = jpos 
            endif 
  170    continue 
         jj = 0 
         len = 0 
!                                                                       
!     eliminate previous rows                                           
!                                                                       
  150    jj = jj+1 
         if (jj  >  lenl) goto 160 
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.         
!-----------------------------------------------------------------------
         jrow = jw(jj) 
         k = jj 
!                                                                       
!     determine smallest column index                                   
!                                                                       
         do 151 j=jj+1,lenl 
            if (jw(j)  <  jrow) then 
               jrow = jw(j) 
               k = j 
            endif 
  151    continue 
!                                                                       
         if (k  /=  jj) then 
!     exchange in jw                                                    
            j = jw(jj) 
            jw(jj) = jw(k) 
            jw(k) = j 
!     exchange in jr                                                    
            jw(n+jrow) = jj 
            jw(n+j) = k 
!     exchange in w                                                     
            s = w(jj) 
            w(jj) = w(k) 
            w(k) = s 
         endif 
!                                                                       
!     zero out element in row by setting resetting jw(n+jrow) to zero.  
!                                                                       
         jw(n+jrow) = 0 
!                                                                       
!     drop term if small                                                
!                                                                       
!         if (abs(w(jj))  <=  tol*tnorm) then                           
!            dropsum = dropsum + w(jj)                                  
!            goto 150                                                   
!         endif                                                         
!                                                                       
!     get the multiplier for row to be eliminated (jrow).               
!                                                                       
         fact = w(jj)*alu(jrow) 
!                                                                       
!     drop term if small                                                
!                                                                       
         if (abs(fact)  <=  tol) then 
            dropsum = dropsum + w(jj) 
            goto 150 
         endif 
!                                                                       
!     combine current row and row jrow                                  
!                                                                       
         do 203 k = ju(jrow), jlu(jrow+1)-1 
            s = fact*alu(k) 
            j = jlu(k) 
            jpos = jw(n+j) 
            if (j  >=  ii) then 
!                                                                       
!     dealing with upper part.                                          
!                                                                       
               if (jpos  ==  0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenu = lenu+1 
                  if (lenu  >  n) goto 995 
                  i = ii+lenu-1 
                  jw(i) = j 
                  jw(n+j) = i 
                  w(i) = - s 
               else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                  w(jpos) = w(jpos) - s 
               endif 
            else 
!                                                                       
!     dealing with lower part.                                          
!                                                                       
               if (jpos  ==  0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenl = lenl+1 
                  if (lenl  >  n) goto 995 
                  jw(lenl) = j 
                  jw(n+j) = lenl 
                  w(lenl) = - s 
               else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                  w(jpos) = w(jpos) - s 
               endif 
            endif 
  203    continue 
         len = len+1 
         w(len) = fact 
         jw(len)  = jrow 
         goto 150 
  160    continue 
!                                                                       
!     reset double-pointer to zero (For U-part only)                    
!                                                                       
         do 308 k=1, lenu 
            jw(n+jw(ii+k-1)) = 0 
  308    continue 
!                                                                       
!     update l-matrix                                                   
!                                                                       
         do 204 k=1, len 
            if (ju0  >  iwk) goto 996 
            alu(ju0) =  w(k) 
            jlu(ju0) =  jw(k) 
            ju0 = ju0+1 
  204    continue 
!                                                                       
!     save pointer to beginning of row ii of U                          
!                                                                       
         ju(ii) = ju0 
!                                                                       
!     go through elements in U-part of w to determine elements to keep  
!                                                                       
         len = 0 
         do k=1, lenu-1 
!            if (abs(w(ii+k))  >  tnorm*tol) then                      
            if (abs(w(ii+k))  >  abs(w(ii))*tol) then 
               len = len+1 
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            else 
               dropsum = dropsum + w(ii+k) 
            endif 
         enddo 
!                                                                       
!     now update u-matrix                                               
!                                                                       
         if (ju0 + len-1  >  iwk) goto 996 
         do 302 k=ii+1,ii+len 
            jlu(ju0) = jw(k) 
            alu(ju0) = w(k) 
            ju0 = ju0+1 
  302    continue 
!                                                                       
!     define diagonal element                                           
!                                                                       
         w(ii) = w(ii) + alph*dropsum 
!                                                                       
!     store inverse of diagonal element of u                            
!                                                                       
         if (w(ii)  ==  0.0) w(ii) = (0.0001 + tol)*tnorm 
!                                                                       
         alu(ii) = 1.0_wp/ w(ii) 
!                                                                       
!     update pointer to beginning of next row of U.                     
!                                                                       
         jlu(ii+1) = ju0 
!-----------------------------------------------------------------------
!     end main loop                                                     
!-----------------------------------------------------------------------
  500 continue 
      ierr = 0 
      return 
!                                                                       
!     incomprehensible error. Matrix must be wrong.                     
!                                                                       
  995 ierr = -1 
      return 
!                                                                       
!     insufficient storage in alu/ jlu arrays for  L / U factors        
!                                                                       
  996 ierr = -2 
      return 
!                                                                       
!     zero row encountered                                              
!                                                                       
  997 ierr = -3 
      return 

      end subroutine ilud                                           

!=============================================================================80

!---------------------------------------------------------------------- 
      subroutine iludp(n,a,ja,ia,alph,droptol,permtol,mbloc,alu,        &
     &     jlu,ju,iwk,w,jw,iperm,ierr)                                  
!-----------------------------------------------------------------------
      implicit none 
      integer n,ja(*),ia(n+1),mbloc,jlu(*),ju(n),jw(2*n),iwk,           &
     &     iperm(2*n),ierr                                              
      real(wp) a(*), alu(*), w(2*n), alph, droptol, permtol 
!----------------------------------------------------------------------*
!                     *** ILUDP preconditioner ***                     *
!    incomplete LU factorization with standard droppoing strategy      *
!    and column pivoting                                               *
!----------------------------------------------------------------------*
! author Yousef Saad -- Aug 1995.                                      *
!----------------------------------------------------------------------*
! on entry:                                                             
!==========                                                             
! n       = integer. The dimension of the matrix A.                     
!                                                                       
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!           ON RETURN THE COLUMNS OF A ARE PERMUTED.                    
!                                                                       
! alph    = diagonal compensation parameter -- the term:                
!                                                                       
!           alph*(sum of all dropped out elements in a given row)       
!                                                                       
!           is added to the diagonal element of U of the factorization  
!           Thus: alph = 0 ---> ~ ILU with threshold,                   
!                 alph = 1 ---> ~ MILU with threshold.                  
!                                                                       
! droptol = tolerance used for dropping elements in L and U.            
!           elements are dropped if they are  <  norm(row) x droptol   
!           row = row being eliminated                                  
!                                                                       
! permtol = tolerance ratio used for determning whether to permute      
!           two columns.  Two columns are permuted only when            
!           abs(a(i,j))*permtol  >  abs(a(i,i))                        
!           [0 --> never permute; good values 0.1 to 0.01]              
!                                                                       
! mbloc   = if desired, permuting can be done only within the diagonal  
!           blocks of size mbloc. Useful for PDE problems with several  
!           degrees of freedom.. If feature not wanted take mbloc=n.    
!                                                                       
! iwk     = integer. The declared lengths of arrays alu and jlu         
!           if iwk is not large enough the code will stop prematurely   
!           with ierr = -2 or ierr = -3 (see below).                    
!                                                                       
! On return:                                                            
!===========                                                            
!                                                                       
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju      = integer array of length n containing the pointers to        
!           the beginning of each row of U in the matrix alu,jlu.       
! iperm   = contains the permutation arrays ..                          
!           iperm(1:n) = old numbers of unknowns                        
!           iperm(n+1:2*n) = reverse permutation = new unknowns.        
!                                                                       
! ierr    = integer. Error message with the following meaning.          
!           ierr  = 0    --> successful return.                         
!           ierr  >  0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.          
!                            (The elimination process has generated a   
!                            row in L or U whose length is  >   n.)    
!           ierr  = -2   --> The L/U matrix overflows the arrays alu,jlu
!           ierr  = -3   --> zero row encountered.                      
!                                                                       
! work arrays:                                                          
!=============                                                          
! jw      = integer work array of length 2*n.                           
! w       = real work array of length 2*n                               
!                                                                       
! Notes:                                                                
! ------                                                                
! IMPORTANT: TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH    
! LU-SOLVE, THE MATRIX A IS PERMUTED ON RETURN. [all column indices are 
! changed]. SIMILARLY FOR THE U MATRIX.                                 
! To permute the matrix back to its original state use the loop:        
!                                                                       
!      do k=ia(1), ia(n+1)-1                                            
!         ja(k) = perm(ja(k))                                           
!      enddo                                                            
!                                                                       
!-----------------------------------------------------------------------
!     local variables                                                   
!                                                                       
      integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,icut 
      real(wp) s,tmp,tnorm,xmax,xmax0,fact,abs,t,dropsum 
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)    
!     and pointer array.                                                
!-----------------------------------------------------------------------
      ju0 = n+2 
      jlu(1) = ju0 
!                                                                       
!  integer double pointer array.                                        
!                                                                       
      do 1 j=1,n 
         jw(n+j)  = 0 
         iperm(j) = j 
         iperm(n+j) = j 
    1 continue 
!-----------------------------------------------------------------------
!     beginning of main loop.                                           
!-----------------------------------------------------------------------
      do 500 ii = 1, n 
         j1 = ia(ii) 
         j2 = ia(ii+1) - 1 
         dropsum = 0.0_wp 
         tnorm = 0.0_wp 
         do 501 k=j1,j2 
            tnorm = tnorm+abs(a(k)) 
  501    continue 
         if (tnorm  ==  0.0) goto 997 
         tnorm = tnorm/(j2-j1+1) 
!                                                                       
!     unpack L-part and U-part of row of A in arrays  w  --             
!                                                                       
         lenu = 1 
         lenl = 0 
         jw(ii) = ii 
         w(ii) = 0.0 
         jw(n+ii) = ii 
!                                                                       
         do 170  j = j1, j2 
            k = iperm(n+ja(j)) 
            t = a(j) 
            if (k  <  ii) then 
               lenl = lenl+1 
               jw(lenl) = k 
               w(lenl) = t 
               jw(n+k) = lenl 
            else if (k  ==  ii) then 
               w(ii) = t 
            else 
               lenu = lenu+1 
               jpos = ii+lenu-1 
               jw(jpos) = k 
               w(jpos) = t 
               jw(n+k) = jpos 
            endif 
  170    continue 
         jj = 0 
         len = 0 
!                                                                       
!     eliminate previous rows                                           
!                                                                       
  150    jj = jj+1 
         if (jj  >  lenl) goto 160 
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.         
!-----------------------------------------------------------------------
         jrow = jw(jj) 
         k = jj 
!                                                                       
!     determine smallest column index                                   
!                                                                       
         do 151 j=jj+1,lenl 
            if (jw(j)  <  jrow) then 
               jrow = jw(j) 
               k = j 
            endif 
  151    continue 
!                                                                       
         if (k  /=  jj) then 
!     exchange in jw                                                    
            j = jw(jj) 
            jw(jj) = jw(k) 
            jw(k) = j 
!     exchange in jr                                                    
            jw(n+jrow) = jj 
            jw(n+j) = k 
!     exchange in w                                                     
            s = w(jj) 
            w(jj) = w(k) 
            w(k) = s 
         endif 
!                                                                       
!     zero out element in row by resetting jw(n+jrow) to zero.          
!                                                                       
         jw(n+jrow) = 0 
!                                                                       
!     drop term if small                                                
!                                                                       
         if (abs(w(jj))  <=  droptol*tnorm) then 
            dropsum = dropsum + w(jj) 
            goto 150 
         endif 
!                                                                       
!     get the multiplier for row to be eliminated: jrow                 
!                                                                       
         fact = w(jj)*alu(jrow) 
!                                                                       
!     combine current row and row jrow                                  
!                                                                       
         do 203 k = ju(jrow), jlu(jrow+1)-1 
            s = fact*alu(k) 
!     new column number                                                 
            j = iperm(n+jlu(k)) 
            jpos = jw(n+j) 
!                                                                       
!     if fill-in element is small then disregard:                       
!                                                                       
            if (j  >=  ii) then 
!                                                                       
!     dealing with upper part.                                          
!                                                                       
               if (jpos  ==  0) then 
!     this is a fill-in element                                         
                  lenu = lenu+1 
                  i = ii+lenu-1 
                  if (lenu  >  n) goto 995 
                  jw(i) = j 
                  jw(n+j) = i 
                  w(i) = - s 
               else 
!     no fill-in element --                                             
                  w(jpos) = w(jpos) - s 
               endif 
            else 
!                                                                       
!     dealing with lower part.                                          
!                                                                       
               if (jpos  ==  0) then 
!     this is a fill-in element                                         
                 lenl = lenl+1 
                 if (lenl  >  n) goto 995 
                 jw(lenl) = j 
                 jw(n+j) = lenl 
                 w(lenl) = - s 
              else 
!     no fill-in element --                                             
                 w(jpos) = w(jpos) - s 
              endif 
           endif 
  203       continue 
        len = len+1 
        w(len) = fact 
        jw(len)  = jrow 
        goto 150 
  160   continue 
!                                                                       
!     reset double-pointer to zero (U-part)                             
!                                                                       
        do 308 k=1, lenu 
           jw(n+jw(ii+k-1)) = 0 
  308       continue 
!                                                                       
!     update L-matrix                                                   
!                                                                       
        do 204 k=1, len 
           if (ju0  >  iwk) goto 996 
           alu(ju0) =  w(k) 
           jlu(ju0) = iperm(jw(k)) 
           ju0 = ju0+1 
  204   continue 
!                                                                       
!     save pointer to beginning of row ii of U                          
!                                                                       
        ju(ii) = ju0 
!                                                                       
!     update u-matrix -- first apply dropping strategy                  
!                                                                       
         len = 0 
         do k=1, lenu-1 
            if (abs(w(ii+k))  >  tnorm*droptol) then 
               len = len+1 
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            else 
               dropsum = dropsum + w(ii+k) 
            endif 
         enddo 
!                                                                       
        imax = ii 
        xmax = abs(w(imax)) 
        xmax0 = xmax 
        icut = ii - 1 + mbloc - mod(ii-1,mbloc) 
!                                                                       
!     determine next pivot --                                           
!                                                                       
        do k=ii+1,ii+len 
           t = abs(w(k)) 
           if (t  >  xmax .and. t*permtol  >  xmax0 .and.             &
     &          jw(k)  <=  icut) then                                   
              imax = k 
              xmax = t 
           endif 
        enddo 
!                                                                       
!     exchange w's                                                      
!                                                                       
        tmp = w(ii) 
        w(ii) = w(imax) 
        w(imax) = tmp 
!                                                                       
!     update iperm and reverse iperm                                    
!                                                                       
        j = jw(imax) 
        i = iperm(ii) 
        iperm(ii) = iperm(j) 
        iperm(j) = i 
!     reverse iperm                                                     
        iperm(n+iperm(ii)) = ii 
        iperm(n+iperm(j)) = j 
!-----------------------------------------------------------------------
        if (len + ju0-1  >  iwk) goto 996 
!                                                                       
!     copy U-part in original coordinates                               
!                                                                       
        do 302 k=ii+1,ii+len 
           jlu(ju0) = iperm(jw(k)) 
           alu(ju0) = w(k) 
           ju0 = ju0+1 
  302       continue 
!                                                                       
!     define diagonal element                                           
!                                                                       
         w(ii) = w(ii) + alph*dropsum 
!                                                                       
!     store inverse of diagonal element of u                            
!                                                                       
        if (w(ii)  ==  0.0) w(ii) = (1.0D-4 + droptol)*tnorm 
!                                                                       
        alu(ii) = 1.0_wp/ w(ii) 
!                                                                       
!     update pointer to beginning of next row of U.                     
!                                                                       
        jlu(ii+1) = ju0 
!-----------------------------------------------------------------------
!     end main loop                                                     
!-----------------------------------------------------------------------
  500 continue 
!                                                                       
!     permute all column indices of LU ...                              
!                                                                       
      do k = jlu(1),jlu(n+1)-1 
         jlu(k) = iperm(n+jlu(k)) 
      enddo 
!                                                                       
!     ...and of A                                                       
!                                                                       
      do k=ia(1), ia(n+1)-1 
         ja(k) = iperm(n+ja(k)) 
      enddo 
!                                                                       
      ierr = 0 
      return 
!                                                                       
!     incomprehensible error. Matrix must be wrong.                     
!                                                                       
  995 ierr = -1 
      return 
!                                                                       
!     insufficient storage in arrays alu, jlu to store factors          
!                                                                       
  996 ierr = -2 
      return 
!                                                                       
!     zero row encountered                                              
!                                                                       
  997 ierr = -3 
      return 

      end subroutine iludp                                           

!=============================================================================80

      subroutine iluk(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr) 

      implicit none 
      integer n 
      real(wp) a(*),alu(*),w(n) 
      integer ja(*),ia(n+1),jlu(*),ju(n),levs(*),jw(3*n),lfil,iwk,ierr 
!----------------------------------------------------------------------*
!     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) *
!----------------------------------------------------------------------*
!                                                                       
! on entry:                                                             
!==========                                                             
! n       = integer. The row dimension of the matrix A. The matrix      
!                                                                       
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!                                                                       
! lfil    = integer. The fill-in parameter. Each element whose          
!           leve-of-fill exceeds lfil during the ILU process is dropped.
!           lfil must be  >=  0                                         
!                                                                       
! iwk     = integer. The minimum length of arrays alu and jlu.
!                                                                       
! On return:                                                            
!===========                                                            
!                                                                       
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju      = integer array of length n containing the pointers to        
!           the beginning of each row of U in the matrix alu,jlu.       
!                                                                       
! levs    = integer (work) array of size iwk -- which contains the      
!           levels of each element in alu, jlu.                         
!                                                                       
! ierr    = integer. Error message with the following meaning.          
!           ierr  = 0    --> successful return.                         
!           ierr  >  0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.          
!                            (The elimination process has generated a   
!                            row in L or U whose length is  >   n.)    
!           ierr  = -2   --> The matrix L overflows the array al.       
!           ierr  = -3   --> The matrix U overflows the array alu.      
!           ierr  = -4   --> Illegal value for lfil.                    
!           ierr  = -5   --> zero row encountered in A or U.            
!                                                                       
! work arrays:                                                          
!=============                                                          
! jw      = integer work array of length 3*n.                           
! w       = real work array of length n                                 
!                                                                       
! Notes/known bugs: This is not implemented efficiently storage-wise.   
!       For example: Only the part of the array levs(*) associated with 
!       the U-matrix is needed in the routine.. So some storage can     
!       be saved if needed. The levels of fills in the LU matrix are    
!       output for information only -- they are not needed by LU-solve. 
!                                                                       
!---------------------------------------------------------------------- 
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]       
! jw(n+1:2n)  stores the nonzero indicator.                             
!                                                                       
! Notes:                                                                
! ------                                                                
! All the diagonal elements of the input matrix must be  nonzero.       
!                                                                       
!----------------------------------------------------------------------*
!     locals                                                            
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,n2,             &
     &     jlev, min                                                    
      real(wp) t, s, fact 
      if (lfil  <  0) goto 998 
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)    
!     and pointer array.                                                
!-----------------------------------------------------------------------
      n2 = n+n 
      ju0 = n+2 
      jlu(1) = ju0 
!                                                                       
!     initialize nonzero indicator array + levs array --                
!                                                                       
      do 1 j=1,2*n 
         jw(j)  = 0 
    1 continue 
!-----------------------------------------------------------------------
!     beginning of main loop.                                           
!-----------------------------------------------------------------------
      do 500 ii = 1, n 
         j1 = ia(ii) 
         j2 = ia(ii+1) - 1 
!                                                                       
!     unpack L-part and U-part of row of A in arrays w                  
!                                                                       
         lenu = 1 
         lenl = 0 
         jw(ii) = ii 
          w(ii) = 0.0 
         jw(n+ii) = ii 
!                                                                       
         do 170  j = j1, j2 
            k = ja(j) 
            t = a(j) 
            if (t  ==  0.0) goto 170 
            if (k  <  ii) then 
               lenl = lenl+1 
               jw(lenl) = k 
               w(lenl) = t 
               jw(n2+lenl) = 0 
               jw(n+k) = lenl 
            else if (k  ==  ii) then 
               w(ii) = t 
               jw(n2+ii) = 0 
            else 
               lenu = lenu+1 
               jpos = ii+lenu-1 
               jw(jpos) = k 
               w(jpos) = t 
               jw(n2+jpos) = 0 
               jw(n+k) = jpos 
            endif 
  170    continue 
!                                                                       
         jj = 0 
!                                                                       
!     eliminate previous rows                                           
!                                                                       
  150    jj = jj+1 
         if (jj  >  lenl) goto 160 
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.         
!-----------------------------------------------------------------------
         jrow = jw(jj) 
         k = jj 
!                                                                       
!     determine smallest column index                                   
!                                                                       
         do 151 j=jj+1,lenl 
            if (jw(j)  <  jrow) then 
               jrow = jw(j) 
               k = j 
            endif 
  151    continue 
!                                                                       
         if (k  /=  jj) then 
!     exchange in jw                                                    
            j = jw(jj) 
            jw(jj) = jw(k) 
            jw(k) = j 
!     exchange in jw(n+  (pointers/ nonzero indicator).                 
            jw(n+jrow) = jj 
            jw(n+j) = k 
!     exchange in jw(n2+  (levels)                                      
            j = jw(n2+jj) 
            jw(n2+jj)  = jw(n2+k) 
            jw(n2+k) = j 
!     exchange in w                                                     
            s = w(jj) 
            w(jj) = w(k) 
            w(k) = s 
         endif 
!                                                                       
!     zero out element in row by resetting jw(n+jrow) to zero.          
!                                                                       
         jw(n+jrow) = 0 
!                                                                       
!     get the multiplier for row to be eliminated (jrow) + its level    
!                                                                       
         fact = w(jj)*alu(jrow) 
         jlev = jw(n2+jj) 
         if (jlev  >  lfil) goto 150 
!                                                                       
!     combine current row and row jrow                                  
!                                                                       
         do 203 k = ju(jrow), jlu(jrow+1)-1 
            s = fact*alu(k) 
            j = jlu(k) 
            jpos = jw(n+j) 
            if (j  >=  ii) then 
!                                                                       
!     dealing with upper part.                                          
!                                                                       
               if (jpos  ==  0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenu = lenu+1 
                  if (lenu  >  n) goto 995 
                  i = ii+lenu-1 
                  jw(i) = j 
                  jw(n+j) = i 
                  w(i) = - s 
                  jw(n2+i) = jlev+levs(k)+1 
               else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                  w(jpos) = w(jpos) - s 
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1) 
               endif 
            else 
!                                                                       
!     dealing with lower part.                                          
!                                                                       
               if (jpos  ==  0) then 
!                                                                       
!     this is a fill-in element                                         
!                                                                       
                  lenl = lenl+1 
                  if (lenl  >  n) goto 995 
                  jw(lenl) = j 
                  jw(n+j) = lenl 
                  w(lenl) = - s 
                  jw(n2+lenl) = jlev+levs(k)+1 
               else 
!                                                                       
!     this is not a fill-in element                                     
!                                                                       
                  w(jpos) = w(jpos) - s 
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1) 
               endif 
            endif 
  203    continue 
         w(jj) = fact 
         jw(jj)  = jrow 
         goto 150 
  160    continue 
!                                                                       
!     reset double-pointer to zero (U-part)                             
!                                                                       
         do 308 k=1, lenu 
            jw(n+jw(ii+k-1)) = 0 
  308    continue 
!                                                                       
!     update l-matrix                                                   
!                                                                       
         do 204 k=1, lenl 
            if (ju0  >  iwk) goto 996 
            if (jw(n2+k)  <=  lfil) then 
               alu(ju0) =  w(k) 
               jlu(ju0) =  jw(k) 
               ju0 = ju0+1 
            endif 
  204    continue 
!                                                                       
!     save pointer to beginning of row ii of U                          
!                                                                       
         ju(ii) = ju0 
!                                                                       
!     update u-matrix                                                   
!                                                                       
         do 302 k=ii+1,ii+lenu-1 
            if (jw(n2+k)  <=  lfil) then 
               jlu(ju0) = jw(k) 
               alu(ju0) = w(k) 
               levs(ju0) = jw(n2+k) 
               ju0 = ju0+1 
            endif 
  302    continue 
                                                                        
         if (w(ii)  ==  0.0) goto 999 
!                                                                       
         alu(ii) = 1.0_wp/ w(ii) 
!                                                                       
!     update pointer to beginning of next row of U.                     
!                                                                       
         jlu(ii+1) = ju0 
!-----------------------------------------------------------------------
!     end main loop                                                     
!-----------------------------------------------------------------------
  500 continue 
      ierr = 0 
      return 
!                                                                       
!     incomprehensible error. Matrix must be wrong.                     
!                                                                       
  995 ierr = -1 
      return 
!                                                                       
!     insufficient storage in L.                                        
!                                                                       
  996 ierr = -2 
      return 
!                                                                       
!     insufficient storage in U.                                        
!                                                                       
  997 ierr = -3 
      return 
!                                                                       
!     illegal lfil entered.                                             
!                                                                       
  998 ierr = -4 
      return 
!                                                                       
!     zero row encountered in A or U.                                   
!                                                                       
  999 ierr = -5 
      write(*,*)'zero row is',ii
      write(*,*)'j_min(ii), j_max(ii),',j1,j2
      do j = j1,j2
        write(*,*)j,a(j)
      enddo
      return 

      end subroutine iluk                                           

!=============================================================================80

        subroutine ilu0(n, a, ja, ia, alu, jlu, ju, iw, ierr) 

        implicit none 

        integer, dimension(*),    intent(in)    ::  ja, ia
        integer, dimension(*),    intent(out)   ::  ju, jlu
        integer, dimension(*),    intent(inout) ::  iw

        integer,                  intent(in)    :: n
        integer,                  intent(out)   :: ierr

        real(wp),  dimension(*), intent(in)  :: a
        real(wp),  dimension(*), intent(out) :: alu

        integer                   :: ju0,js,j,jcol,jf,jm,i,ii
        integer                   :: jj,jw,jrow
        real(wp)               :: tl


!------------------ right preconditioner ------------------------------*
!                    ***   ilu(0) preconditioner.   ***                *
!----------------------------------------------------------------------*
! Note that this has been coded in such a way that it can be used       
! with pgmres. Normally, since the data structure of the L+U matrix is  
! the same as that the A matrix, savings can be made. In fact with      
! some definitions (not correct for general sparse matrices) all we     
! need in addition to a, ja, ia is an additional diagonal.              
! ILU0 is not recommended for serious problems. It is only provided     
! here for comparison purposes.                                         
!-----------------------------------------------------------------------
!                                                                       
! on entry:                                                             
!---------                                                              
! n       = dimension of matrix                                         
! a, ja,                                                                
! ia      = original matrix in compressed sparse row storage.           
!                                                                       
! on return:                                                            
!-----------                                                            
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju          = pointer to the diagonal elements in alu, jlu.           
!                                                                       
! ierr          = integer indicating error code on return               
!             ierr = 0 --> normal return                                
!             ierr = k --> code encountered a zero pivot at step k.     
! work arrays:                                                          
!-------------                                                          
! iw            = integer work array of length n.                       
!------------                                                           
! IMPORTANT                                                             
!-----------                                                            
! it is assumed that the the elements in the input matrix are stored    
!    in such a way that in each row the lower part comes first and      
!    then the upper part. To get the correct ILU factorization, it is   
!    also necessary to have the elements of L sorted by increasing      
!    column number. It may therefore be necessary to sort the           
!    elements of a, ja, ia prior to calling ilu0. This can be           
!    achieved by transposing the matrix twice using csrcsc.             
!                                                                       
!-----------------------------------------------------------------------
        ju0 = n+2 
        jlu(1) = ju0 
!                                                                       
! initialize work vector to zero's                                      
!                                                                       
        do 31 i=1, n 
           iw(i) = 0 
   31   continue 
!                                                                       
! main loop                                                             
!                                                                       
        do 500 ii = 1, n 
           js = ju0 
!                                                                       
! generating row number ii of L and U.                                  
!                                                                       
           do 100 j=ia(ii),ia(ii+1)-1 
!                                                                       
!     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.    
!                                                                       
              jcol = ja(j) 
              if (jcol  ==  ii) then 
                 alu(ii) = a(j) 
                 iw(jcol) = ii 
                 ju(ii)  = ju0 
              else 
                 alu(ju0) = a(j) 
                 jlu(ju0) = ja(j) 
                 iw(jcol) = ju0 
                 ju0 = ju0+1 
              endif 
  100      continue 
           jlu(ii+1) = ju0 
           jf = ju0-1 
           jm = ju(ii)-1 
!                                                                       
!     exit if diagonal element is reached.                              
!                                                                       
           do 150 j=js, jm 
              jrow = jlu(j) 
              tl = alu(j)*alu(jrow) 
              alu(j) = tl 
!                                                                       
!     perform  linear combination                                       
!                                                                       
              do 140 jj = ju(jrow), jlu(jrow+1)-1 
                 jw = iw(jlu(jj)) 
                 if (jw  /=  0) alu(jw) = alu(jw) - tl*alu(jj) 
  140         continue 
  150      continue 
!                                                                       
!     invert  and store diagonal element.                               
!                                                                       
           if (alu(ii)  ==  0.0_wp) goto 600 
           alu(ii) = 1.0_wp/alu(ii) 
!                                                                       
!     reset pointer iw to zero                                          
!                                                                       
           iw(ii) = 0 
           do 201 i = js, jf 
  201         iw(jlu(i)) = 0 
  500      continue 
           ierr = 0 
           return 
!                                                                       
!     zero pivot :                                                      
!                                                                       
  600      ierr = ii 
!                                                                       
           return 

      end subroutine ilu0                                           

!=============================================================================80

        subroutine milu0(n, a, ja, ia, alu, jlu, ju, iw, ierr) 

        implicit none 

        integer, dimension(*),    intent(in)    ::  ja, ia
        integer, dimension(*),    intent(out)   ::  ju, jlu
        integer, dimension(*),    intent(inout) ::  iw

        integer,                  intent(in)    :: n
        integer,                  intent(out)   :: ierr

        real(wp),  dimension(*), intent(in)  :: a
        real(wp),  dimension(*), intent(out) :: alu

        integer                   :: ju0,js,j,jcol,jf,jm,i,ii
        integer                   :: jj,jw,jrow
        real(wp)               :: tl, s

!       here is the original dimensions for this subroutine.  The ones above
!       were added based on the ilu0 routine (very similar). 
!       THIS ROUTINE HAS NOT BEEN VERIFIED!

!       implicit real(wp) (a-h,o-z) 
!       real(wp) a(*), alu(*) 
!       integer ja(*), ia(*), ju(*), jlu(*), iw(*) 

!----------------------------------------------------------------------*
!                *** simple milu(0) preconditioner. ***                *
!----------------------------------------------------------------------*
! Note that this has been coded in such a way that it can be used       
! with pgmres. Normally, since the data structure of a, ja, ia is       
! the same as that of a, ja, ia, savings can be made. In fact with      
! some definitions (not correct for general sparse matrices) all we     
! need in addition to a, ja, ia is an additional diagonal.              
! Ilu0 is not recommended for serious problems. It is only provided     
! here for comparison purposes.                                         
!-----------------------------------------------------------------------
!                                                                       
! on entry:                                                             
!----------                                                             
! n       = dimension of matrix                                         
! a, ja,                                                                
! ia      = original matrix in compressed sparse row storage.           
!                                                                       
! on return:                                                            
!----------                                                             
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in       
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix 
!           contains the i-th row of L (excluding the diagonal entry=1) 
!           followed by the i-th row of U.                              
!                                                                       
! ju          = pointer to the diagonal elements in alu, jlu.           
!                                                                       
! ierr          = integer indicating error code on return               
!             ierr = 0 --> normal return                                
!             ierr = k --> code encountered a zero pivot at step k.     
! work arrays:                                                          
!-------------                                                          
! iw            = integer work array of length n.                       
!------------                                                           
! Note (IMPORTANT):                                                     
!-----------                                                            
! it is assumed that the the elements in the input matrix are ordered   
!    in such a way that in each row the lower part comes first and      
!    then the upper part. To get the correct ILU factorization, it is   
!    also necessary to have the elements of L ordered by increasing     
!    column number. It may therefore be necessary to sort the           
!    elements of a, ja, ia prior to calling milu0. This can be          
!    achieved by transposing the matrix twice using csrcsc.             
!-----------------------------------------------------------            
          ju0 = n+2 
          jlu(1) = ju0 
! initialize work vector to zero's                                      
        do 31 i=1, n 
   31         iw(i) = 0 
!                                                                       
!-------------- MAIN LOOP ----------------------------------            
!                                                                       
        do 500 ii = 1, n 
           js = ju0 
!                                                                       
! generating row number ii or L and U.                                  
!                                                                       
           do 100 j=ia(ii),ia(ii+1)-1 
!                                                                       
!     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.    
!                                                                       
              jcol = ja(j) 
              if (jcol  ==  ii) then 
                 alu(ii) = a(j) 
                 iw(jcol) = ii 
                 ju(ii)  = ju0 
              else 
                 alu(ju0) = a(j) 
                 jlu(ju0) = ja(j) 
                 iw(jcol) = ju0 
                 ju0 = ju0+1 
              endif 
  100      continue 
           jlu(ii+1) = ju0 
           jf = ju0-1 
           jm = ju(ii)-1 
!     s accumulates fill-in values                                      
           s = 0.0_wp 
           do 150 j=js, jm 
              jrow = jlu(j) 
              tl = alu(j)*alu(jrow) 
              alu(j) = tl 
!-----------------------perform linear combination --------             
              do 140 jj = ju(jrow), jlu(jrow+1)-1 
                 jw = iw(jlu(jj)) 
                 if (jw  /=  0) then 
                       alu(jw) = alu(jw) - tl*alu(jj) 
                    else 
                       s = s + tl*alu(jj) 
                    endif 
  140         continue 
  150      continue 
!----------------------- invert and store diagonal element.             
           alu(ii) = alu(ii)-s 
           if (alu(ii)  ==  0.0_wp) goto 600 
           alu(ii) = 1.0_wp/alu(ii) 
!----------------------- reset pointer iw to zero                       
           iw(ii) = 0 
           do 201 i = js, jf 
  201         iw(jlu(i)) = 0 
  500      continue 
           ierr = 0 
           return 
!     zero pivot :                                                      
  600      ierr = ii 
           return 

      end subroutine milu0                                           

!=============================================================================80

        subroutine ilu0_csr(n, nnz, a, ja, ia, iau, aI, jaI, iaI, iauI, iw, ierr)

        implicit none 

        integer,                  intent(in)    :: n, nnz
        integer,                  intent(out)   :: ierr

        integer, dimension(*),    intent(in)    ::  ja, ia, iau
        integer, dimension(*),    intent(out)   ::  jaI, iaI, iauI

        real(wp),  dimension(*), intent(in)  :: a
        real(wp),  dimension(*), intent(out) :: aI

        integer, dimension(*),    intent(inout) ::  iw


        integer                   :: i,k
        integer                   :: j,jj,j1,j2,jw,jrow
        real(wp)               :: t


!------------------ right preconditioner ------------------------------*
!                    ***   ilu(0) preconditioner.   ***                *
!----------------------------------------------------------------------*
! Note that this has been coded in such a way that it can be used       
! with pgmres. Normally, since the data structure of the L+U matrix is  
! the same as that the A matrix, savings can be made. In fact with      
! some definitions (not correct for general sparse matrices) all we     
! need in addition to a, ja, ia is an additional diagonal.              
! ILU0 is not recommended for serious problems. It is only provided     
! here for comparison purposes.                                         
!-----------------------------------------------------------------------
!                                                                       
! on entry:                                                             
!---------                                                              
! n           = dimension of matrix                                         
! a           = original matrix in compressed sparse row storag
! ja, ia, iau 
!                                                                       
! on return:                                                            
!-----------                                                            
!                                                                       
! ierr          = integer indicating error code on return               
!             ierr = 0 --> normal return                                
!             ierr = k --> code encountered a zero pivot at step k.     
! work arrays:                                                          
!-------------                                                          
! iw            = integer work array of length n.                       
!------------                                                           
! IMPORTANT                                                             
!-----------                                                            
! it is assumed that the the elements in the input matrix are stored    
!    in such a way that in each row the lower part comes first and      
!    then the upper part. To get the correct ILU factorization, it is   
!    also necessary to have the elements of L sorted by increasing      
!    column number. It may therefore be necessary to sort the           
!    elements of a, ja, ia prior to calling ilu0. This can be           
!    achieved by transposing the matrix twice using csrcsc.             
!                                                                       
!-----------------------------------------------------------------------

!                                                                       
! copy original matrix "a" into initial L U matrix s                                      
!                                                                       
        do i=1, nnz 
           aI(i) = a(i)
        enddo

        do i=1, nnz 
           jaI(i) = ja(i)
        enddo
        do i=1, n+1 
           iaI(i)  = ia(i)
        enddo
        do i=1, n
           iauI(i) = iau(i)
        enddo
!                                                                       
! initialize work vector to zero's                                      
!                                                                       
        do i=1, n 
           iw(i) = 0 
        enddo

!                                                                       
!-------------- MAIN LOOP ----------------------------------
                                                                        
        do k = 1, n 

!------------------------ k = row number -------------------

          j1 = ia(k)
          j2 = ia(k+1)-1
          do j = j1, j2
            iw(ja(j)) = j
          enddo

          j=j1
150       jrow = ja(j)

!----------------------- Exit if diagonal element is reached.

          if (jrow >= k) goto 200
   
!----------------------- Compute the multiplier for jrow.

          t      = aI(j) * aI(iau(jrow))

          aI(j) = t
!                                                                       
!     perform  linear combination                                       
!                                                                       
          do jj = iau(jrow)+1, ia(jrow+1)-1 
             jw = iw(ja(jj)) 
             if (jw  /=  0) aI(jw) = aI(jw) - t*aI(jj) 
          enddo

          j = j + 1
          if (j <= j2) goto 150

  200     continue

!                                                                       
!     invert  and store diagonal element.                               
!                                                                       

          if (jrow /= k ) goto 600 

          aI(j) = 1.0_wp/aI(j) 
!                                                                       
!     reset pointer iw to zero                                          
!                                                                       
          do i = j1, j2 
            iw(ja(i)) = 0 
          enddo

        enddo

        ierr = 0 
        return 
!                                                                       
!   zero pivot :                                                      
!                                                                       
  600   ierr = k 
!                                                                       
        return 

      end subroutine ilu0_csr

!=============================================================================80

       subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout,        &
     &  aa, ja, ia, alu, jlu, ju, ierr)

       use eispack_module,   only : epslon
       use blas_module,      only : daxpy,ddot,dnrm2

       implicit none

       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
       real(wp) vv(n,*), rhs(n), sol(n), aa(*), alu(*)

       integer                   :: n1,its,j,ro,i,i1,k,k1,ii,jj,kmax
       real(wp),          intent(in)        :: eps
       real(wp)                             :: eps1,gam

!!!-----------------------------------------------------------------------
!!       subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout,        &
!!     &                    aa, ja, ia, alu, jlu, ju, ierr)               
!!!-----------------------------------------------------------------------
!!       implicit real(wp) (a-h,o-z) 
!!       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n) 
!!       real(wp) vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps 
!----------------------------------------------------------------------*
!                                                                      *
!                 *** ILUT - Preconditioned GMRES ***                  *
!                                                                      *
!----------------------------------------------------------------------*
! This is a simple version of the ILUT preconditioned GMRES algorithm. *
! The ILUT preconditioner uses a dual strategy for dropping elements   *
! instead  of the usual level of-fill-in approach. See details in ILUT *
! subroutine documentation. PGMRES uses the L and U matrices generated *
! from the subroutine ILUT to precondition the GMRES algorithm.        *
! The preconditioning is applied to the right. The stopping criterion  *
! utilized is based simply on reducing the residual norm by epsilon.   *
! This preconditioning is more reliable than ilu0 but requires more    *
! storage. It seems to be much less prone to difficulties related to   *
! strong nonsymmetries in the matrix. We recommend using a nonzero tol *
! (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
! lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
! more reliable the code is. Efficiency may also be much improved.     *
! Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
! Gaussian elimination without pivoting.                               *
!                                                                      *
! ILU(0) and MILU(0) are also provided for comparison purposes         *
! USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
! then call pgmres.                                                    *
!----------------------------------------------------------------------*
! Coded by Y. Saad - This version dated May, 7, 1990.                  *
!----------------------------------------------------------------------*
! parameters                                                           *
!-----------                                                           *
! on entry:                                                            *
!==========                                                            *
!                                                                      *
! n     == integer. The dimension of the matrix.                       *
! im    == size of krylov subspace:  should not exceed 50 in this      *
!          version (can be reset by changing parameter command for     *
!          kmax below)                                                 *
! rhs   == real vector of length n containing the right hand side.     *
!          Destroyed on return.                                        *
! sol   == real vector of length n containing an initial guess to the  *
!          solution on input. approximate solution on output           *
! eps   == tolerance for stopping criterion. process is stopped        *
!          as soon as ( ||.|| is the euclidean norm):                  *
!          || current residual||/||initial residual|| <= eps           *
! maxits== maximum number of iterations allowed                        *
! iout  == output unit number number for printing intermediate results *
!          if (iout  <=  0) nothing is printed out.                    *
!                                                                      *
! aa, ja,                                                              *
! ia    == the input matrix in compressed sparse row format:           *
!          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
!          ja(1:nnz) = corresponding column indices.                   *
!          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
!          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
!                                                                      *
! alu,jlu== A matrix stored in Modified Sparse Row format containing   *
!           the L and U factors, as computed by subroutine ilut.       *
!                                                                      *
! ju     == integer array of length n containing the pointers to       *
!           the beginning of each row of U in alu, jlu as computed     *
!           by subroutine ILUT.                                        *
!                                                                      *
! on return:                                                           *
!==========                                                            *
! sol   == contains an approximate solution (upon successful return).  *
! ierr  == integer. Error message with the following meaning.          *
!          ierr = 0 --> successful return.                             *
!          ierr = 1 --> convergence not achieved in itmax iterations.  *
!          ierr =-1 --> the initial guess seems to be the exact        *
!                       solution (initial residual computed was zero)  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! work arrays:                                                         *
!=============                                                         *
! vv    == work array of length  n x (im+1) (used to store the Arnoli  *
!          basis)                                                      *
!----------------------------------------------------------------------*
! subroutines called :                                                 *
! amux   : SPARSKIT routine to do the matrix by vector multiplication  *
!          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
! lusol : combined forward and backward solves (Preconditioning ope.) * 
! BLAS1  routines.                                                     *
!----------------------------------------------------------------------*
       parameter (kmax=50) 
       real(wp) hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t 
!-------------------------------------------------------------          
! arnoldi size should not exceed kmax=50 in this version..              
! to reset modify paramter kmax accordingly.                            
!-------------------------------------------------------------          
!      data epsmac/1.d-16/ 
       n1 = n + 1 
       its = 0 
!-------------------------------------------------------------          
! outer loop starts here..                                              
!-------------- compute initial residual vector --------------          
       call amux (n, sol, vv, aa, ja, ia) 
       do 21 j=1,n 
          vv(j,1) = rhs(j) - vv(j,1) 
   21  continue 
!-------------------------------------------------------------          
   20  ro = dnrm2(n, vv, 1) 
       if (iout  >  0 .and. its  ==  0)                                &
     &      write(iout, 199) its, ro                                    
       if (ro  ==  0.0_wp) goto 999 
       t = 1.0_wp/ ro 
       do 210 j=1, n 
          vv(j,1) = vv(j,1)*t 
  210  continue 
       if (its  ==  0) eps1=eps*ro 
!     ** initialize 1-st term  of rhs of hessenberg system..            
       rs(1) = ro 
       i = 0 
    4  i=i+1 
       its = its + 1 
       i1 = i + 1 
       call lusol (n, vv(1,i), rhs, alu, jlu, ju) 
       call amux (n, rhs, vv(1,i1), aa, ja, ia) 
!-----------------------------------------                              
!     modified gram - schmidt...                                        
!-----------------------------------------                              
       do 55 j=1, i 
          t = ddot(n, vv(1,j),1,vv(1,i1),1) 
          hh(j,i) = t 
          call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1) 
   55  continue 
       t = dnrm2(n, vv(1,i1), 1) 
       hh(i1,i) = t 
       if ( t  ==  0.0_wp) goto 58 
       t = 1.0_wp/t 
       do 57  k=1,n 
          vv(k,i1) = vv(k,i1)*t 
   57  continue 
!                                                                       
!     done with modified gram schimd and arnoldi step..                 
!     now  update factorization of hh                                   
!                                                                       
   58  if (i  ==  1) goto 121 
!--------perfrom previous transformations  on i-th column of h          
       do 66 k=2,i 
          k1 = k-1 
          t = hh(k1,i) 
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i) 
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i) 
   66  continue 
  121  gam = sqrt(hh(i,i)**2 + hh(i1,i)**2) 
!                                                                       
!     if gamma is zero then any small value will do...                  
!     will affect only residual estimate                                
!                                                                       
       if (gam  ==  0.0_wp) gam = epslon(1.0_wp)
!                                                                       
!     get  next plane rotation                                          
!                                                                       
       c(i) = hh(i,i)/gam 
       s(i) = hh(i1,i)/gam 
       rs(i1) = -s(i)*rs(i) 
       rs(i) =  c(i)*rs(i) 
!                                                                       
!     detrermine residual norm and test for convergence-                
!                                                                       
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i) 
       ro = abs(rs(i1)) 
  131  format(1h ,2e14.4) 
       if (iout  >  0)                                                 &
     &      write(iout, 199) its, ro                                    
       if (i  <  im .and. (ro  >  eps1))  goto 4 
!                                                                       
!     now compute solution. first solve upper triangular system.        
!                                                                       
       rs(i) = rs(i)/hh(i,i) 
       do 30 ii=2,i 
          k=i-ii+1 
          k1 = k+1 
          t=rs(k) 
          do 40 j=k1,i 
             t = t-hh(k,j)*rs(j) 
   40     continue 
          rs(k) = t/hh(k,k) 
   30  continue 
!                                                                       
!     form linear combination of v(*,i)'s to get solution               
!                                                                       
       t = rs(1) 
       do 15 k=1, n 
          rhs(k) = vv(k,1)*t 
   15  continue 
       do 16 j=2, i 
          t = rs(j) 
          do 161 k=1, n 
             rhs(k) = rhs(k)+t*vv(k,j) 
  161     continue 
   16  continue 
!                                                                       
!     call preconditioner.                                              
!                                                                       
       call lusol (n, rhs, rhs, alu, jlu, ju) 
       do 17 k=1, n 
          sol(k) = sol(k) + rhs(k) 
   17  continue 
!                                                                       
!     restart outer loop  when necessary                                
!                                                                       
       if (ro  <=  eps1) goto 990 
       if (its  >=  maxits) goto 991 
!                                                                       
!     else compute residual vector and continue..                       
!                                                                       
       do 24 j=1,i 
          jj = i1-j+1 
          rs(jj-1) = -s(jj-1)*rs(jj) 
          rs(jj) = c(jj-1)*rs(jj) 
   24  continue 
       do 25  j=1,i1 
          t = rs(j) 
          if (j  ==  1)  t = t-1.0_wp 
          call daxpy (n, t, vv(1,j), 1,  vv, 1) 
   25  continue 
  199  format('   its =', i4, ' res. norm =', d20.6) 
!     restart outer loop.                                               
       goto 20 
  990  ierr = 0 
       return 
  991  ierr = 1 
       return 
  999  continue 
       ierr = -1 
       return 

      end subroutine pgmres                                           

!=============================================================================80

        subroutine lusol(n, y, x, alu, jlu, ju) 

        implicit none 

        integer,      dimension(*), intent(in)    :: jlu, ju
        integer,                    intent(in)    :: n

        real(wp),  dimension(*), intent(in)    :: alu, y
        real(wp),  dimension(*), intent(out)   :: x

!-----------------------------------------------------------------------
!                                                                       
! This routine solves the system (LU) x = y,                            
! given an LU decomposition of a matrix stored in (alu, jlu, ju)        
! MODIFIED sparse row format                                            
!                                                                       
!-----------------------------------------------------------------------
! on entry:                                                             
! n   = dimension of system                                             
! y   = the right-hand-side vector                                      
! alu, jlu, ju                                                          
!     = the LU matrix as provided from the ILU routines.                
!                                                                       
! on return                                                             
! x   = solution of LU x = y.                                           
!-----------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)         
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
        integer i,k 
                                                                        
!       forward solve                                                         
                                                                        
        do i = 1, n 
           x(i) = y(i) 
           do k=jlu(i),ju(i)-1 
              x(i) = x(i) - alu(k)* x(jlu(k)) 
           enddo
        enddo

!       backward solve.                                                   
                                                                        
        do i = n, 1, -1 
           do k=ju(i),jlu(i+1)-1 
              x(i) = x(i) - alu(k)*x(jlu(k)) 
           enddo
           x(i) = alu(i)*x(i) 
        enddo
!                                                                       
        return 

        end subroutine lusol                                           

!=============================================================================80

        subroutine lusol_M_rhs(n, m, y, x, alu, jlu, ju) 

        implicit none 

        integer,                      intent(in)  :: n, m
        integer,      dimension(  *), intent(in)  :: jlu, ju

        real(wp),  dimension(  *), intent(in)  :: alu
        real(wp),  dimension(m,n), intent(in)  :: y

        real(wp),  dimension(m,n), intent(out) :: x

!-----------------------------------------------------------------------
!                                                                       
! This routine solves the system (LU) x = y,                            
! given an LU decomposition of a matrix stored in (alu, jlu, ju)        
! MODIFIED sparse row format                                            
!                                                                       
!-----------------------------------------------------------------------
! on entry:                                                             
! n   = dimension of system                                             
! y   = the right-hand-side vector                                      
! alu, jlu, ju                                                          
!     = the LU matrix as provided from the ILU routines.                
!                                                                       
! on return                                                             
! x   = solution of LU x = y.                                           
!-----------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)         
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
        integer i,k 
                                                                        
!       forward solve                                                         
                                                                        
        do i = 1, n 
           x(:,i) = y(:,i) 
           do k=jlu(i),ju(i)-1 
              x(:,i) = x(:,i) - alu(k)* x(:,jlu(k)) 
           enddo
        enddo

!       backward solve.                                                   
                                                                        
        do i = n, 1, -1 
           do k=ju(i),jlu(i+1)-1 
              x(:,i) = x(:,i) - alu(k)*x(:,jlu(k)) 
           enddo
           x(:,i) = alu(i)*x(:,i) 
        enddo
!                                                                       
        return 

        end subroutine lusol_M_rhs

!=============================================================================80

        subroutine B_lusol(n, neq, nnz, y, x, alu, jlu, ju) 

        implicit none 

        integer,      dimension(*), intent(in)    :: jlu, ju
        integer,                    intent(in)    :: n, neq, nnz

        real(wp),  dimension(neq,neq,nnz), intent(in)    :: alu 
        real(wp),  dimension(neq,nnz)    , intent(in)    :: y 
        real(wp),  dimension(neq,nnz)    , intent(out)   :: x 

!-----------------------------------------------------------------------
!                                                                       
! This routine solves the system (LU) x = y,                            
! given an LU decomposition of a matrix stored in (alu, jlu, ju)        
! MODIFIED sparse row format                                            
!                                                                       
!-----------------------------------------------------------------------
! on entry:                                                             
! n   = dimension of system                                             
! y   = the right-hand-side vector                                      
! alu, jlu, ju                                                          
!     = the LU matrix as provided from the ILU routines.                
!                                                                       
! on return                                                             
! x   = solution of LU x = y.                                           
!-----------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)         
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
        integer i,k , L1, L2
        real(wp),  dimension(neq)    :: t 
                                                                        
!       forward solve                                                         
                                                                        
        do i = 1, n 

          do L1 = 1,neq
            x(L1,i) = y(L1,i) 
            do k=jlu(i),ju(i)-1 

              do L2 = 1,neq
                x(L1,i) = x(L1,i) - alu(L1,L2,k)* x(L2,jlu(k)) 
              enddo

            enddo
          enddo

        enddo

!       backward solve.                                                   
                                                                        
        do i = n, 1, -1 

          do L1 = 1,neq
            do k=ju(i),jlu(i+1)-1 

              do L2 = 1,neq
                x(L1,i) = x(L1,i) - alu(L1,L2,k)*x(L2,jlu(k)) 
              enddo

            enddo
          enddo

          do L1 = 1,neq
            t(L1) = 0.0_wp
            do L2 = 1,neq
              t(L1) = t(L1) + alu(L1,L2,i)*x(L2,i) 
            enddo
          enddo
          do L1 = 1,neq
            x(L1,i) = t(L1)
          enddo

        enddo
!                                                                       
        return 

        end subroutine B_lusol

!=============================================================================80

        subroutine lusol_csr(n, y, x, a, ja, ia, iau) 

        implicit none 

        integer,      dimension(*), intent(in)    :: ja, ia, iau
        integer,                    intent(in)    :: n

        real(wp),  dimension(*), intent(in)    :: a, y
        real(wp),  dimension(*), intent(out)   :: x

!-----------------------------------------------------------------------
!                                                                       
!-----------------------------------------------------------------------
!   solves    L x = y ; L = lower unit triang.       /  CSR format
!   followed by 
!             U x = x    U = unit upper triangular.  /  CSR format
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! a,     =  values of matrix in CSR
! ja,    =  columns of matrix for each "a"
! ia,    =  pointer to first element of each row in compressed sparse row
! iau,   =  pointer to diagonal element
!
! On return:
!-----------
!       x  = The solution of  L U x  = y.
!--------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lusol (n, x, x, a, ja, ia)         
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
                                                                        
      integer        :: i, k

!                                                                       
! forward solve                                                         
!                                                                       
      do i = 1, n
         x(i) = y(i)
         do k = ia(i), iau(i)-1
            x(i) = x(i)-a(k)*x(ja(k))
         enddo
      enddo

!                                                                       
!     backward solve.                                                   
!                                                                       
      do i = n,1,-1
         do k = iau(i)+1, ia(i+1)-1
            x(i) = x(i) - a(k)*x(ja(k))
         enddo
         x(i) = x(i) * a(iau(i))
      enddo

      end subroutine lusol_csr

!=============================================================================80

        subroutine lutsol(n, y, x, alu, jlu, ju) 
        real(wp) x(n), y(n), alu(*) 
        integer n, jlu(*), ju(*) 
!-----------------------------------------------------------------------
!                                                                       
! This routine solves the system  Transp(LU) x = y,                     
! given an LU decomposition of a matrix stored in (alu, jlu, ju)        
! modified sparse row format. Transp(M) is the transpose of M.          
!-----------------------------------------------------------------------
! on entry:                                                             
! n   = dimension of system                                             
! y   = the right-hand-side vector                                      
! alu, jlu, ju                                                          
!     = the LU matrix as provided from the ILU routines.                
!                                                                       
! on return                                                             
! x   = solution of transp(LU) x = y.                                   
!-----------------------------------------------------------------------
!                                                                       
! Note: routine is in place: call lutsol (n, x, x, alu, jlu, ju)        
!       will solve the system with rhs x and overwrite the result on x .
!                                                                       
!-----------------------------------------------------------------------
! local variables                                                       
!                                                                       
        integer i,k 
!                                                                       
        do 10 i = 1, n 
           x(i) = y(i) 
   10   continue 
!                                                                       
! forward solve (with U^T)                                              
!                                                                       
        do 20 i = 1, n 
           x(i) = x(i) * alu(i) 
           do 30 k=ju(i),jlu(i+1)-1 
              x(jlu(k)) = x(jlu(k)) - alu(k)* x(i) 
   30      continue 
   20   continue 
!                                                                       
!     backward solve (with L^T)                                         
!                                                                       
        do 40 i = n, 1, -1 
           do 50 k=jlu(i),ju(i)-1 
              x(jlu(k)) = x(jlu(k)) - alu(k)*x(i) 
   50      continue 
   40   continue 
!                                                                       
          return 

      end subroutine lutsol                                           

!=============================================================================80

!-----------------------------------------------------------------------
        subroutine qsplit(a,ind,n,ncut) 

        real(wp) a(n) 
        integer ind(n), n, ncut 
        integer mid,j
!-----------------------------------------------------------------------
!     does a quick-sort split of a real array.                          
!     on input a(1:n). is a real array                                  
!     on output a(1:n) is permuted such that its elements satisfy:      
!                                                                       
!     abs(a(i))  >=  abs(a(ncut)) for i  <  ncut and                   
!     abs(a(i))  <=  abs(a(ncut)) for i  >  ncut                       
!                                                                       
!     ind(1:n) is an integer array which permuted in the same way as a(*
!-----------------------------------------------------------------------
        real(wp) tmp, abskey 
        integer itmp, first, last 
!-----                                                                  
        first = 1 
        last = n 
        if (ncut  <  first .or. ncut  >  last) return 
!                                                                       
!     outer loop -- while mid  /=  ncut do                              
!                                                                       
    1   mid = first 
        abskey = abs(a(mid)) 
        do 2 j=first+1, last 
           if (abs(a(j))  >  abskey) then 
              mid = mid+1 
!     interchange                                                       
              tmp = a(mid) 
              itmp = ind(mid) 
              a(mid) = a(j) 
              ind(mid) = ind(j) 
              a(j)  = tmp 
              ind(j) = itmp 
           endif 
    2   continue 
!                                                                       
!     interchange                                                       
!                                                                       
        tmp = a(mid) 
        a(mid) = a(first) 
        a(first)  = tmp 
!                                                                       
        itmp = ind(mid) 
        ind(mid) = ind(first) 
        ind(first) = itmp 
!                                                                       
!     test for while loop                                               
!                                                                       
        if (mid  ==  ncut) return 
        if (mid  >  ncut) then 
           last = mid-1 
        else 
           first = mid+1 
        endif 
        goto 1 

      end subroutine qsplit                                           

!=============================================================================80

      subroutine msrcsr (n,a,ja,ao,jao,iao,iauo,wk,iwk) 

      integer,   intent(in)  :: n

      integer ja(*),jao(*),iao(n+1),iauo(n),iwk(n+1) 
      real(wp) a(*),ao(*),wk(n) 
!-----------------------------------------------------------------------
!       Modified - Sparse Row  to   Compressed Sparse Row               
!                                                                       
!-----------------------------------------------------------------------
! converts a compressed matrix using a separated diagonal               
! (modified sparse row format) in the Compressed Sparse Row             
! format.                                                               
! does not check for zero elements in the diagonal.                     
!                                                                       
!                                                                       
! on entry :                                                            
!---------                                                              
! n          = row dimension of matrix                                  
! a, ja      = sparse matrix in msr sparse storage format               
!              see routine csrmsr for details on data structure         
!                                                                       
! on return :                                                           
!-----------                                                            
!                                                                       
! ao,jao,iao = output matrix in csr format.                             
!                                                                       
! work arrays:                                                          
!------------                                                           
! wk       = real work array of length n                                
! iwk      = integer work array of length n+1                           
!                                                                       
! notes:                                                                
!   The original version of this was NOT in place, but has              
!   been modified by adding the vector iwk to be in place.              
!   The original version had ja instead of iwk everywhere in            
!   loop 500.  Modified  Sun 29 May 1994 by R. Bramley (Indiana).       
!                                                                       
!-----------------------------------------------------------------------

      integer     ::  i,ii,j,k,iptr,idiag
      logical added 

      do 1 i=1,n 
         wk(i) = a(i) 
         iwk(i) = ja(i) 
    1 continue 
      iwk(n+1) = ja(n+1) 
      iao(1) = 1 
      iptr = 1 
!---------                                                              
      do 500 ii=1,n 
         added = .false. 
         idiag = iptr + (iwk(ii+1)-iwk(ii)) 
         do 100 k=iwk(ii),iwk(ii+1)-1 
            j = ja(k) 
            if (j  <  ii) then 
               ao(iptr) = a(k) 
               jao(iptr) = j 
               iptr = iptr+1 
            elseif (added) then 
               ao(iptr) = a(k) 
               jao(iptr) = j 
               iptr = iptr+1 
            else 
! add diag element - only reserve a position for it.                    
               idiag = iptr 
               iptr = iptr+1 
               added = .true. 
!     then other element                                                
               ao(iptr) = a(k) 
               jao(iptr) = j 
               iptr = iptr+1 
            endif 
  100    continue 
         ao(idiag) = wk(ii) 
         jao(idiag) = ii 
         if (.not. added) iptr = iptr+1 
         iao(ii+1) = iptr 
         do k=iao(ii),iao(ii+1)-1
           if(jao(k) == ii) iauo(ii) = k
         enddo
  500 continue 
      return 

      end subroutine msrcsr                                           

      end module ilut_module
