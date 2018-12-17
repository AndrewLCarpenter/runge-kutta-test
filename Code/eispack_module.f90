      module eispack_module

      use precision_vars      , only : wp

      implicit none; save

      private

      public ::  rgg , rg, rs, svd, hqr, qsortd, minfit
      public ::  elmhes, eltran, epslon
      public ::  Eigen_Sym

!     integer, parameter               :: dp = 8    ! double precision

      character (len=*), parameter :: eispack_module_scm_id =   &
      '$Id: eispack_module.f90 19469 2006-06-09 11:32:26Z kleb $'

      contains


!     balanc, balbak, cdiv, elmhes, eltran,
!     hqr1, hqr2, hqr, imtqlv, minfit,
!     tqlrat, qzhes, qzit, qzval, qzvec, qsortd,
!     rebak, reduc, rg, rgg, rs,
!     rsg, rsm, svd, tinvit, tql1,
!     tql2, trbak1, tred1, tred2


!=============================================================================80

      subroutine balanc(nm,n,a,low,igh,scale) 
!                                                                       
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc 
      real(wp) a(nm,n),scale(n) 
      real(wp) c,f,g,r,s,b2,radix 
      logical noconv 
!                                                                       
!     this subroutine is a translation of the algol procedure balance,  
!     num. math. 13, 293-304(1969) by parlett and reinsch.              
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).   
!                                                                       
!     this subroutine balances a real matrix and isolates               
!     eigenvalues whenever possible.                                    
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        a contains the input matrix to be balanced.                    
!                                                                       
!     on output                                                         
!                                                                       
!        a contains the balanced matrix.                                
!                                                                       
!        low and igh are two integers such that a(i,j)                  
!          is equal to zero if                                          
!           (1) i is greater than j and                                 
!           (2) j=1,...,low-1 or i=igh+1,...,n.                         
!                                                                       
!        scale contains information determining the                     
!           permutations and scaling factors used.                      
!                                                                       
!     suppose that the principal submatrix in rows low through igh      
!     has been balanced, that p(j) denotes the index interchanged       
!     with j during the permutation step, and that the elements         
!     of the diagonal matrix used are denoted by d(i,j).  then          
!        scale(j) = p(j),    for j = 1,...,low-1                        
!                 = d(j,j),      j = low,...,igh                        
!                 = p(j)         j = igh+1,...,n.                       
!     the order in which the interchanges are made is n to igh+1,       
!     then 1 to low-1.                                                  
!                                                                       
!     note that 1 is returned for igh if igh is zero formally.          
!                                                                       
!     the algol procedure exc contained in balance appears in           
!     balanc  in line.  (note that the algol roles of identifiers       
!     k,l have been reversed.)                                          
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      radix = 16.0_wp 
!                                                                       
      b2 = radix * radix 
      k = 1 
      l = n 
      go to 100 
!     .......... in-line procedure for row and                          
!                column exchange ..........                             
   20 scale(m) = j 
      if (j  ==  m) go to 50 
!                                                                       
      do 30 i = 1, l 
         f = a(i,j) 
         a(i,j) = a(i,m) 
         a(i,m) = f 
   30 continue 
!                                                                       
      do 40 i = k, n 
         f = a(j,i) 
         a(j,i) = a(m,i) 
         a(m,i) = f 
   40 continue 
!                                                                       
   50 go to (80,130), iexc 
!     .......... search for rows isolating an eigenvalue                
!                and push them down ..........                          
   80 if (l  ==  1) go to 280 
      l = l - 1 
!     .......... for j=l step -1 until 1 do -- ..........               
  100 do 120 jj = 1, l 
         j = l + 1 - jj 
!                                                                       
         do 110 i = 1, l 
            if (i  ==  j) go to 110 
            if (a(j,i)  /=  0.0_wp) go to 120 
  110    continue 
!                                                                       
         m = l 
         iexc = 1 
         go to 20 
  120 continue 
!                                                                       
      go to 140 
!     .......... search for columns isolating an eigenvalue             
!                and push them left ..........                          
  130 k = k + 1 
!                                                                       
  140 do 170 j = k, l 
!                                                                       
         do 150 i = k, l 
            if (i  ==  j) go to 150 
            if (a(i,j)  /=  0.0_wp) go to 170 
  150    continue 
!                                                                       
         m = k 
         iexc = 2 
         go to 20 
  170 continue 
!     .......... now balance the submatrix in rows k to l ..........    
      do 180 i = k, l 
  180 scale(i) = 1.0_wp 
!     .......... iterative loop for norm reduction ..........           
  190 noconv = .false. 
!                                                                       
      do 270 i = k, l 
         c = 0.0_wp 
         r = 0.0_wp 
!                                                                       
         do 200 j = k, l 
            if (j  ==  i) go to 200 
            c = c + abs(a(j,i)) 
            r = r + abs(a(i,j)) 
  200    continue 
!     .......... guard against zero c or r due to underflow ..........  
         if (c  ==  0.0_wp .or. r  ==  0.0_wp) go to 270 
         g = r / radix 
         f = 1.0_wp 
         s = c + r 
  210    if (c  >=  g) go to 220 
         f = f * radix 
         c = c * b2 
         go to 210 
  220    g = r * radix 
  230    if (c  <  g) go to 240 
         f = f / radix 
         c = c / b2 
         go to 230 
!     .......... now balance ..........                                 
  240    if ((c + r) / f  >=  0.95_wp * s) go to 270 
         g = 1.0_wp / f 
         scale(i) = scale(i) * f 
         noconv = .true. 
!                                                                       
         do 250 j = k, n 
  250    a(i,j) = a(i,j) * g 
!                                                                       
         do 260 j = 1, l 
  260    a(j,i) = a(j,i) * f 
!                                                                       
  270 continue 
!                                                                       
      if (noconv) go to 190 
!                                                                       
  280 low = k 
      igh = l 
      return 
      end subroutine balanc

!=============================================================================80

      subroutine balbak(nm,n,low,igh,scale,m,z) 
!                                                                       
      integer i,j,k,m,n,ii,nm,igh,low 
      real(wp) scale(n),z(nm,m) 
      real(wp) s 
!                                                                       
!     this subroutine is a translation of the algol procedure balbak,   
!     num. math. 13, 293-304(1969) by parlett and reinsch.              
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).   
!                                                                       
!     this subroutine forms the eigenvectors of a real general          
!     matrix by back transforming those of the corresponding            
!     balanced matrix determined by  balanc.                            
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        low and igh are integers determined by  balanc.                
!                                                                       
!        scale contains information determining the permutations        
!          and scaling factors used by  balanc.                         
!                                                                       
!        m is the number of columns of z to be back transformed.        
!                                                                       
!        z contains the real and imaginary parts of the eigen-          
!          vectors to be back transformed in its first m columns.       
!                                                                       
!     on output                                                         
!                                                                       
!        z contains the real and imaginary parts of the                 
!          transformed eigenvectors in its first m columns.             
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (m  ==  0) go to 200 
      if (igh  ==  low) go to 120 
!                                                                       
      do 110 i = low, igh 
         s = scale(i) 
!     .......... left hand eigenvectors are back transformed            
!                if the foregoing statement is replaced by              
!                s=1.0_wp/scale(i). ..........                           
         do 100 j = 1, m 
  100    z(i,j) = z(i,j) * s 
!                                                                       
  110 continue 
!     ......... for i=low-1 step -1 until 1,                            
!               igh+1 step 1 until n do -- ..........                   
  120 do 140 ii = 1, n 
         i = ii 
         if (i  >=  low .and. i  <=  igh) go to 140 
         if (i  <  low) i = low - ii 
         k = scale(i) 
         if (k  ==  i) go to 140 
!                                                                       
         do 130 j = 1, m 
            s = z(i,j) 
            z(i,j) = z(k,j) 
            z(k,j) = s 
  130    continue 
!                                                                       
  140 continue 
!                                                                       
  200 return 
      end subroutine balbak

!=============================================================================80

      subroutine cdiv(ar,ai,br,bi,cr,ci) 
      real(wp) ar,ai,br,bi,cr,ci 
!                                                                       
!     complex division, (cr,ci) = (ar,ai)/(br,bi)                       
!                                                                       
      real(wp) s,ars,ais,brs,bis 
      s = abs(br) + abs(bi) 
      ars = ar/s 
      ais = ai/s 
      brs = br/s 
      bis = bi/s 
      s = brs**2 + bis**2 
      cr = (ars*brs + ais*bis)/s 
      ci = (ais*brs - ars*bis)/s 
      return 
      end subroutine cdiv

!=============================================================================80

      subroutine elmhes(nm,n,low,igh,a,int) 
!                                                                       
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1 
      real(wp) a(nm,n) 
      real(wp) x,y 
      integer int(igh) 
!                                                                       
!     this subroutine is a translation of the algol procedure elmhes,   
!     num. math. 12, 349-368(1968) by martin and wilkinson.             
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).   
!                                                                       
!     given a real general matrix, this subroutine                      
!     reduces a submatrix situated in rows and columns                  
!     low through igh to upper hessenberg form by                       
!     stabilized elementary similarity transformations.                 
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        low and igh are integers determined by the balancing           
!          subroutine  balanc.  if  balanc  has not been used,          
!          set low=1, igh=n.                                            
!                                                                       
!        a contains the input matrix.                                   
!                                                                       
!     on output                                                         
!                                                                       
!        a contains the hessenberg matrix.  the multipliers             
!          which were used in the reduction are stored in the           
!          remaining triangle under the hessenberg matrix.              
!                                                                       
!        int contains information on the rows and columns               
!          interchanged in the reduction.                               
!          only elements low through igh are used.                      
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      la = igh - 1 
      kp1 = low + 1 
      if (la  <  kp1) go to 200 
!                                                                       
      do 180 m = kp1, la 
         mm1 = m - 1 
         x = 0.0_wp 
         i = m 
!                                                                       
         do 100 j = m, igh 
            if (abs(a(j,mm1))  <=  abs(x)) go to 100 
            x = a(j,mm1) 
            i = j 
  100    continue 
!                                                                       
         int(m) = i 
         if (i  ==  m) go to 130 
!     .......... interchange rows and columns of a ..........           
         do 110 j = mm1, n 
            y = a(i,j) 
            a(i,j) = a(m,j) 
            a(m,j) = y 
  110    continue 
!                                                                       
         do 120 j = 1, igh 
            y = a(j,i) 
            a(j,i) = a(j,m) 
            a(j,m) = y 
  120    continue 
!     .......... end interchange ..........                             
  130    if (x  ==  0.0_wp) go to 180 
         mp1 = m + 1 
!                                                                       
         do 160 i = mp1, igh 
            y = a(i,mm1) 
            if (y  ==  0.0_wp) go to 160 
            y = y / x 
            a(i,mm1) = y 
!                                                                       
            do 140 j = m, n 
  140       a(i,j) = a(i,j) - y * a(m,j) 
!                                                                       
            do 150 j = 1, igh 
  150       a(j,m) = a(j,m) + y * a(j,i) 
!                                                                       
  160    continue 
!                                                                       
  180 continue 
!                                                                       
  200 return 
      end subroutine elmhes

!=============================================================================80

      subroutine eltran(nm,n,low,igh,a,int,z) 
!                                                                       
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1 
      real(wp) a(nm,igh),z(nm,n) 
      integer int(igh) 
!                                                                       
!     this subroutine is a translation of the algol procedure elmtrans, 
!     num. math. 16, 181-204(1970) by peters and wilkinson.             
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).   
!                                                                       
!     this subroutine accumulates the stabilized elementary             
!     similarity transformations used in the reduction of a             
!     real general matrix to upper hessenberg form by  elmhes.          
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        low and igh are integers determined by the balancing           
!          subroutine  balanc.  if  balanc  has not been used,          
!          set low=1, igh=n.                                            
!                                                                       
!        a contains the multipliers which were used in the              
!          reduction by  elmhes  in its lower triangle                  
!          below the subdiagonal.                                       
!                                                                       
!        int contains information on the rows and columns               
!          interchanged in the reduction by  elmhes.                    
!          only elements low through igh are used.                      
!                                                                       
!     on output                                                         
!                                                                       
!        z contains the transformation matrix produced in the           
!          reduction by  elmhes.                                        
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
!     .......... initialize z to identity matrix ..........             
      do 80 j = 1, n 
!                                                                       
         do 60 i = 1, n 
   60    z(i,j) = 0.0_wp 
!                                                                       
         z(j,j) = 1.0_wp 
   80 continue 
!                                                                       
      kl = igh - low - 1 
      if (kl  <  1) go to 200 
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........      
      do 140 mm = 1, kl 
         mp = igh - mm 
         mp1 = mp + 1 
!                                                                       
         do 100 i = mp1, igh 
  100    z(i,mp) = a(i,mp-1) 
!                                                                       
         i = int(mp) 
         if (i  ==  mp) go to 140 
!                                                                       
         do 130 j = mp, igh 
            z(mp,j) = z(i,j) 
            z(i,j) = 0.0_wp 
  130    continue 
!                                                                       
         z(i,mp) = 1.0_wp 
  140 continue 
!                                                                       
  200 return 
      end subroutine eltran

!=============================================================================80

      real(wp) function epslon (x) 
      real(wp) x 
!                                                                       
!     estimate unit roundoff in quantities of size x.                   
!                                                                       
      real(wp) a,b,c,eps 
!                                                                       
!     this program should function properly on all systems              
!     satisfying the following two assumptions,                         
!        1.  the base used in representing floating point               
!            numbers is not a power of three.                           
!        2.  the quantity  a  in statement 10 is represented to         
!            the accuracy used in floating point variables              
!            that are stored in memory.                                 
!     the statement number 10 and the go to 10 are intended to          
!     force optimizing compilers to generate code satisfying            
!     assumption 2.                                                     
!     under these assumptions, it should be true that,                  
!            a  is not exactly equal to four-thirds,                    
!            b  has a zero for its last bit or digit,                   
!            c  is not exactly equal to one,                            
!            eps  measures the separation of 1.0 from                   
!                 the next larger floating point number.                
!     the developers of eispack would appreciate being informed         
!     about any systems where these assumptions do not hold.            
!                                                                       
!     this version dated 4/6/83.                                        
!                                                                       
      a = 4.0_wp/3.0_wp 
   10 b = a - 1.0_wp 
      c = b + b + b 
      eps = abs(c-1.0_wp) 
      if (eps  ==  0.0_wp) go to 10 
      epslon = eps*abs(x) 
      return 
      end function epslon

!=============================================================================80

      subroutine hqr1(nm,n,low,igh,h,wr,wi,ierr) 
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)   
!                                                                       
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr 
      real(wp) h(nm,n),wr(n),wi(n) 
      real(wp) p,q,r,s,t,w,x,y,zz,norm,tst1,tst2 
      logical notlas 
!                                                                       
!     this subroutine is a translation of the algol procedure hqr,      
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.    
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).   
!                                                                       
!     this subroutine finds the eigenvalues of a real                   
!     upper hessenberg matrix by the qr method.                         
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        low and igh are integers determined by the balancing           
!          subroutine  balanc.  if  balanc  has not been used,          
!          set low=1, igh=n.                                            
!                                                                       
!        h contains the upper hessenberg matrix.  information about     
!          the transformations used in the reduction to hessenberg      
!          form by  elmhes  or  orthes, if performed, is stored         
!          in the remaining triangle under the hessenberg matrix.       
!                                                                       
!     on output                                                         
!                                                                       
!        h has been destroyed.  therefore, it must be saved             
!          before calling  hqr  if subsequent calculation and           
!          back transformation of eigenvectors is to be performed.      
!                                                                       
!        wr and wi contain the real and imaginary parts,                
!          respectively, of the eigenvalues.  the eigenvalues           
!          are unordered except that complex conjugate pairs            
!          of values appear consecutively with the eigenvalue           
!          having the positive imaginary part first.  if an             
!          error exit is made, the eigenvalues should be correct        
!          for indices ierr+1,...,n.                                    
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the limit of 30*n iterations is exhausted      
!                     while the j-th eigenvalue is being sought.        
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated september 1989.                                
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      norm = 0.0_wp 
      k = 1 
!     .......... store roots isolated by balanc                         
!                and compute matrix norm ..........                     
      do 50 i = 1, n 
!                                                                       
         do 40 j = k, n 
   40    norm = norm + abs(h(i,j)) 
!                                                                       
         k = i 
         if (i  >=  low .and. i  <=  igh) go to 50 
         wr(i) = h(i,i) 
         wi(i) = 0.0_wp 
   50 continue 
!                                                                       
      en = igh 
      t = 0.0_wp 
      itn = 30*n 
!     .......... search for next eigenvalues ..........                 
   60 if (en  <  low) go to 1001 
      its = 0 
      na = en - 1 
      enm2 = na - 1 
!     .......... look for single small sub-diagonal element             
!                for l=en step -1 until low do -- ..........            
   70 do 80 ll = low, en 
         l = en + low - ll 
         if (l  ==  low) go to 100 
         s = abs(h(l-1,l-1)) + abs(h(l,l)) 
         if (s  ==  0.0_wp) s = norm 
         tst1 = s 
         tst2 = tst1 + abs(h(l,l-1)) 
         if (tst2  ==  tst1) go to 100 
   80 continue 
!     .......... form shift ..........                                  
  100 x = h(en,en) 
      if (l  ==  en) go to 270 
      y = h(na,na) 
      w = h(en,na) * h(na,en) 
      if (l  ==  na) go to 280 
      if (itn  ==  0) go to 1000 
      if (its  /=  10 .and. its  /=  20) go to 130 
!     .......... form exceptional shift ..........                      
      t = t + x 
!                                                                       
      do 120 i = low, en 
  120 h(i,i) = h(i,i) - x 
!                                                                       
      s = abs(h(en,na)) + abs(h(na,enm2)) 
      x = 0.75_wp * s 
      y = x 
      w = -0.4375_wp * s * s 
  130 its = its + 1 
      itn = itn - 1 
!     .......... look for two consecutive small                         
!                sub-diagonal elements.                                 
!                for m=en-2 step -1 until l do -- ..........            
      do 140 mm = l, enm2 
         m = enm2 + l - mm 
         zz = h(m,m) 
         r = x - zz 
         s = y - zz 
         p = (r * s - w) / h(m+1,m) + h(m,m+1) 
         q = h(m+1,m+1) - zz - r - s 
         r = h(m+2,m+1) 
         s = abs(p) + abs(q) + abs(r) 
         p = p / s 
         q = q / s 
         r = r / s 
         if (m  ==  l) go to 150 
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1))) 
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r)) 
         if (tst2  ==  tst1) go to 150 
  140 continue 
!                                                                       
  150 mp2 = m + 2 
!                                                                       
      do 160 i = mp2, en 
         h(i,i-2) = 0.0_wp 
         if (i  ==  mp2) go to 160 
         h(i,i-3) = 0.0_wp 
  160 continue 
!     .......... double qr step involving rows l to en and              
!                columns m to en ..........                             
      do 260 k = m, na 
         notlas = k  /=  na 
         if (k  ==  m) go to 170 
         p = h(k,k-1) 
         q = h(k+1,k-1) 
         r = 0.0_wp 
         if (notlas) r = h(k+2,k-1) 
         x = abs(p) + abs(q) + abs(r) 
         if (x  ==  0.0_wp) go to 260 
         p = p / x 
         q = q / x 
         r = r / x 
  170    s = sign(sqrt(p*p+q*q+r*r),p) 
         if (k  ==  m) go to 180 
         h(k,k-1) = -s * x 
         go to 190 
  180    if (l  /=  m) h(k,k-1) = -h(k,k-1) 
  190    p = p + s 
         x = p / s 
         y = q / s 
         zz = r / s 
         q = q / p 
         r = r / p 
         if (notlas) go to 225 
!     .......... row modification ..........                            
         do 200 j = k, EN 
            p = h(k,j) + q * h(k+1,j) 
            h(k,j) = h(k,j) - p * x 
            h(k+1,j) = h(k+1,j) - p * y 
  200    continue 
!                                                                       
         j = min0(en,k+3) 
!     .......... column modification ..........                         
         do 210 i = L, j 
            p = x * h(i,k) + y * h(i,k+1) 
            h(i,k) = h(i,k) - p 
            h(i,k+1) = h(i,k+1) - p * q 
  210    continue 
         go to 255 
  225    continue 
!     .......... row modification ..........                            
         do 230 j = k, EN 
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j) 
            h(k,j) = h(k,j) - p * x 
            h(k+1,j) = h(k+1,j) - p * y 
            h(k+2,j) = h(k+2,j) - p * zz 
  230    continue 
!                                                                       
         j = min0(en,k+3) 
!     .......... column modification ..........                         
         do 240 i = L, j 
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2) 
            h(i,k) = h(i,k) - p 
            h(i,k+1) = h(i,k+1) - p * q 
            h(i,k+2) = h(i,k+2) - p * r 
  240    continue 
  255    continue 
!                                                                       
  260 continue 
!                                                                       
      go to 70 
!     .......... one root found ..........                              
  270 wr(en) = x + t 
      wi(en) = 0.0_wp 
      en = na 
      go to 60 
!     .......... two roots found ..........                             
  280 p = (y - x) / 2.0_wp 
      q = p * p + w 
      zz = sqrt(abs(q)) 
      x = x + t 
      if (q  <  0.0_wp) go to 320 
!     .......... real pair ..........                                   
      zz = p + sign(zz,p) 
      wr(na) = x + zz 
      wr(en) = wr(na) 
      if (zz  /=  0.0_wp) wr(en) = x - w / zz 
      wi(na) = 0.0_wp 
      wi(en) = 0.0_wp 
      go to 330 
!     .......... complex pair ..........                                
  320 wr(na) = x + p 
      wr(en) = x + p 
      wi(na) = zz 
      wi(en) = -zz 
  330 en = enm2 
      go to 60 
!     .......... set error -- all eigenvalues have not                  
!                converged after 30*n iterations ..........             
 1000 ierr = en 
 1001 return 
      end subroutine hqr1

!=============================================================================80

      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr) 
!                                                                       
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,                      &
     &        igh,itn,its,low,mp2,enm2,ierr                             
      real(wp) h(nm,n),wr(n),wi(n),z(nm,n) 
      real(wp) p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2 
      logical notlas 
!                                                                       
!     this subroutine is a translation of the algol procedure hqr2,     
!     num. math. 16, 181-204(1970) by peters and wilkinson.             
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).   
!                                                                       
!     this subroutine finds the eigenvalues and eigenvectors            
!     of a real upper hessenberg matrix by the qr method.  the          
!     eigenvectors of a real general matrix can also be found           
!     if  elmhes  and  eltran  or  orthes  and  ortran  have            
!     been used to reduce this general matrix to hessenberg form        
!     and to accumulate the similarity transformations.                 
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        low and igh are integers determined by the balancing           
!          subroutine  balanc.  if  balanc  has not been used,          
!          set low=1, igh=n.                                            
!                                                                       
!        h contains the upper hessenberg matrix.                        
!                                                                       
!        z contains the transformation matrix produced by  eltran       
!          after the reduction by  elmhes, or by  ortran  after the     
!          reduction by  orthes, if performed.  if the eigenvectors     
!          of the hessenberg matrix are desired, z must contain the     
!          identity matrix.                                             
!                                                                       
!     on output                                                         
!                                                                       
!        h has been destroyed.                                          
!                                                                       
!        wr and wi contain the real and imaginary parts,                
!          respectively, of the eigenvalues.  the eigenvalues           
!          are unordered except that complex conjugate pairs            
!          of values appear consecutively with the eigenvalue           
!          having the positive imaginary part first.  if an             
!          error exit is made, the eigenvalues should be correct        
!          for indices ierr+1,...,n.                                    
!                                                                       
!        z contains the real and imaginary parts of the eigenvectors.   
!          if the i-th eigenvalue is real, the i-th column of z         
!          contains its eigenvector.  if the i-th eigenvalue is complex 
!          with positive imaginary part, the i-th and (i+1)-th          
!          columns of z contain the real and imaginary parts of its     
!          eigenvector.  the eigenvectors are unnormalized.  if an      
!          error exit is made, none of the eigenvectors has been found. 
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the limit of 30*n iterations is exhausted      
!                     while the j-th eigenvalue is being sought.        
!                                                                       
!     calls cdiv for complex division.                                  
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      norm = 0.0_wp 
      k = 1 
!     .......... store roots isolated by balanc                         
!                and compute matrix norm ..........                     
      do 50 i = 1, n 
!                                                                       
         do 40 j = k, n 
   40    norm = norm + abs(h(i,j)) 
!                                                                       
         k = i 
         if (i  >=  low .and. i  <=  igh) go to 50 
         wr(i) = h(i,i) 
         wi(i) = 0.0_wp 
   50 continue 
!                                                                       
      en = igh 
      t = 0.0_wp 
      itn = 30*n 
!     .......... search for next eigenvalues ..........                 
   60 if (en  <  low) go to 340 
      its = 0 
      na = en - 1 
      enm2 = na - 1 
!     .......... look for single small sub-diagonal element             
!                for l=en step -1 until low do -- ..........            
   70 do 80 ll = low, en 
         l = en + low - ll 
         if (l  ==  low) go to 100 
         s = abs(h(l-1,l-1)) + abs(h(l,l)) 
         if (s  ==  0.0_wp) s = norm 
         tst1 = s 
         tst2 = tst1 + abs(h(l,l-1)) 
         if (tst2  ==  tst1) go to 100 
   80 continue 
!     .......... form shift ..........                                  
  100 x = h(en,en) 
      if (l  ==  en) go to 270 
      y = h(na,na) 
      w = h(en,na) * h(na,en) 
      if (l  ==  na) go to 280 
      if (itn  ==  0) go to 1000 
      if (its  /=  10 .and. its  /=  20) go to 130 
!     .......... form exceptional shift ..........                      
      t = t + x 
!                                                                       
      do 120 i = low, en 
  120 h(i,i) = h(i,i) - x 
!                                                                       
      s = abs(h(en,na)) + abs(h(na,enm2)) 
      x = 0.75_wp * s 
      y = x 
      w = -0.4375_wp * s * s 
  130 its = its + 1 
      itn = itn - 1 
!     .......... look for two consecutive small                         
!                sub-diagonal elements.                                 
!                for m=en-2 step -1 until l do -- ..........            
      do 140 mm = l, enm2 
         m = enm2 + l - mm 
         zz = h(m,m) 
         r = x - zz 
         s = y - zz 
         p = (r * s - w) / h(m+1,m) + h(m,m+1) 
         q = h(m+1,m+1) - zz - r - s 
         r = h(m+2,m+1) 
         s = abs(p) + abs(q) + abs(r) 
         p = p / s 
         q = q / s 
         r = r / s 
         if (m  ==  l) go to 150 
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1))) 
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r)) 
         if (tst2  ==  tst1) go to 150 
  140 continue 
!                                                                       
  150 mp2 = m + 2 
!                                                                       
      do 160 i = mp2, en 
         h(i,i-2) = 0.0_wp 
         if (i  ==  mp2) go to 160 
         h(i,i-3) = 0.0_wp 
  160 continue 
!     .......... double qr step involving rows l to en and              
!                columns m to en ..........                             
      do 260 k = m, na 
         notlas = k  /=  na 
         if (k  ==  m) go to 170 
         p = h(k,k-1) 
         q = h(k+1,k-1) 
         r = 0.0_wp 
         if (notlas) r = h(k+2,k-1) 
         x = abs(p) + abs(q) + abs(r) 
         if (x  ==  0.0_wp) go to 260 
         p = p / x 
         q = q / x 
         r = r / x 
  170    s = sign(sqrt(p*p+q*q+r*r),p) 
         if (k  ==  m) go to 180 
         h(k,k-1) = -s * x 
         go to 190 
  180    if (l  /=  m) h(k,k-1) = -h(k,k-1) 
  190    p = p + s 
         x = p / s 
         y = q / s 
         zz = r / s 
         q = q / p 
         r = r / p 
         if (notlas) go to 225 
!     .......... row modification ..........                            
         do 200 j = k, n 
            p = h(k,j) + q * h(k+1,j) 
            h(k,j) = h(k,j) - p * x 
            h(k+1,j) = h(k+1,j) - p * y 
  200    continue 
!                                                                       
         j = min0(en,k+3) 
!     .......... column modification ..........                         
         do 210 i = 1, j 
            p = x * h(i,k) + y * h(i,k+1) 
            h(i,k) = h(i,k) - p 
            h(i,k+1) = h(i,k+1) - p * q 
  210    continue 
!     .......... accumulate transformations ..........                  
         do 220 i = low, igh 
            p = x * z(i,k) + y * z(i,k+1) 
            z(i,k) = z(i,k) - p 
            z(i,k+1) = z(i,k+1) - p * q 
  220    continue 
         go to 255 
  225    continue 
!     .......... row modification ..........                            
         do 230 j = k, n 
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j) 
            h(k,j) = h(k,j) - p * x 
            h(k+1,j) = h(k+1,j) - p * y 
            h(k+2,j) = h(k+2,j) - p * zz 
  230    continue 
!                                                                       
         j = min0(en,k+3) 
!     .......... column modification ..........                         
         do 240 i = 1, j 
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2) 
            h(i,k) = h(i,k) - p 
            h(i,k+1) = h(i,k+1) - p * q 
            h(i,k+2) = h(i,k+2) - p * r 
  240    continue 
!     .......... accumulate transformations ..........                  
         do 250 i = low, igh 
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2) 
            z(i,k) = z(i,k) - p 
            z(i,k+1) = z(i,k+1) - p * q 
            z(i,k+2) = z(i,k+2) - p * r 
  250    continue 
  255    continue 
!                                                                       
  260 continue 
!                                                                       
      go to 70 
!     .......... one root found ..........                              
  270 h(en,en) = x + t 
      wr(en) = h(en,en) 
      wi(en) = 0.0_wp 
      en = na 
      go to 60 
!     .......... two roots found ..........                             
  280 p = (y - x) / 2.0_wp 
      q = p * p + w 
      zz = sqrt(abs(q)) 
      h(en,en) = x + t 
      x = h(en,en) 
      h(na,na) = y + t 
      if (q  <  0.0_wp) go to 320 
!     .......... real pair ..........                                   
      zz = p + sign(zz,p) 
      wr(na) = x + zz 
      wr(en) = wr(na) 
      if (zz  /=  0.0_wp) wr(en) = x - w / zz 
      wi(na) = 0.0_wp 
      wi(en) = 0.0_wp 
      x = h(en,na) 
      s = abs(x) + abs(zz) 
      p = x / s 
      q = zz / s 
      r = sqrt(p*p+q*q) 
      p = p / r 
      q = q / r 
!     .......... row modification ..........                            
      do 290 j = na, n 
         zz = h(na,j) 
         h(na,j) = q * zz + p * h(en,j) 
         h(en,j) = q * h(en,j) - p * zz 
  290 continue 
!     .......... column modification ..........                         
      do 300 i = 1, en 
         zz = h(i,na) 
         h(i,na) = q * zz + p * h(i,en) 
         h(i,en) = q * h(i,en) - p * zz 
  300 continue 
!     .......... accumulate transformations ..........                  
      do 310 i = low, igh 
         zz = z(i,na) 
         z(i,na) = q * zz + p * z(i,en) 
         z(i,en) = q * z(i,en) - p * zz 
  310 continue 
!                                                                       
      go to 330 
!     .......... complex pair ..........                                
  320 wr(na) = x + p 
      wr(en) = x + p 
      wi(na) = zz 
      wi(en) = -zz 
  330 en = enm2 
      go to 60 
!     .......... all roots found.  backsubstitute to find               
!                vectors of upper triangular form ..........            
  340 if (norm  ==  0.0_wp) go to 1001 
!     .......... for en=n step -1 until 1 do -- ..........              
      do 800 nn = 1, n 
         en = n + 1 - nn 
         p = wr(en) 
         q = wi(en) 
         na = en - 1 
         if (q) 710, 600, 800 
!     .......... real vector ..........                                 
  600    m = en 
         h(en,en) = 1.0_wp 
         if (na  ==  0) go to 800 
!     .......... for i=en-1 step -1 until 1 do -- ..........            
         do 700 ii = 1, na 
            i = en - ii 
            w = h(i,i) - p 
            r = 0.0_wp 
!                                                                       
            do 610 j = m, en 
  610       r = r + h(i,j) * h(j,en) 
!                                                                       
            if (wi(i)  >=  0.0_wp) go to 630 
            zz = w 
            s = r 
            go to 700 
  630       m = i 
            if (wi(i)  /=  0.0_wp) go to 640 
            t = w 
            if (t  /=  0.0_wp) go to 635 
               tst1 = norm 
               t = tst1 
  632          t = 0.01_wp * t 
               tst2 = norm + t 
               if (tst2  >  tst1) go to 632 
  635       h(i,en) = -r / t 
            go to 680 
!     .......... solve real equations ..........                        
  640       x = h(i,i+1) 
            y = h(i+1,i) 
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) 
            t = (x * s - zz * r) / q 
            h(i,en) = t 
            if (abs(x)  <=  abs(zz)) go to 650 
            h(i+1,en) = (-r - w * t) / x 
            go to 680 
  650       h(i+1,en) = (-s - y * t) / zz 
!                                                                       
!     .......... overflow control ..........                            
  680       t = abs(h(i,en)) 
            if (t  ==  0.0_wp) go to 700 
            tst1 = t 
            tst2 = tst1 + 1.0_wp/tst1 
            if (tst2  >  tst1) go to 700 
            do 690 j = i, en 
               h(j,en) = h(j,en)/t 
  690       continue 
!                                                                       
  700    continue 
!     .......... end real vector ..........                             
         go to 800 
!     .......... complex vector ..........                              
  710    m = na 
!     .......... last vector component chosen imaginary so that         
!                eigenvector matrix is triangular ..........            
         if (abs(h(en,na))  <=  abs(h(na,en))) go to 720 
         h(na,na) = q / h(en,na) 
         h(na,en) = -(h(en,en) - p) / h(en,na) 
         go to 730 
  720    call cdiv(0.0_wp,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en)) 
  730    h(en,na) = 0.0_wp 
         h(en,en) = 1.0_wp 
         enm2 = na - 1 
         if (enm2  ==  0) go to 800 
!     .......... for i=en-2 step -1 until 1 do -- ..........            
         do 795 ii = 1, enm2 
            i = na - ii 
            w = h(i,i) - p 
            ra = 0.0_wp 
            sa = 0.0_wp 
!                                                                       
            do 760 j = m, en 
               ra = ra + h(i,j) * h(j,na) 
               sa = sa + h(i,j) * h(j,en) 
  760       continue 
!                                                                       
            if (wi(i)  >=  0.0_wp) go to 770 
            zz = w 
            r = ra 
            s = sa 
            go to 795 
  770       m = i 
            if (wi(i)  /=  0.0_wp) go to 780 
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en)) 
            go to 790 
!     .......... solve complex equations ..........                     
  780       x = h(i,i+1) 
            y = h(i+1,i) 
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q 
            vi = (wr(i) - p) * 2.0_wp * q 
            if (vr  /=  0.0_wp .or. vi  /=  0.0_wp) go to 784 
               tst1 = norm * (abs(w) + abs(q) + abs(x)               &
     &                      + abs(y) + abs(zz))                       
               vr = tst1 
  783          vr = 0.01_wp * vr 
               tst2 = tst1 + vr 
               if (tst2  >  tst1) go to 783 
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,              &
     &                h(i,na),h(i,en))                                  
            if (abs(x)  <=  abs(zz) + abs(q)) go to 785 
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x 
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x 
            go to 790 
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,                   &
     &                h(i+1,na),h(i+1,en))                              
!                                                                       
!     .......... overflow control ..........                            
  790       t = max(abs(h(i,na)), abs(h(i,en))) 
            if (t  ==  0.0_wp) go to 795 
            tst1 = t 
            tst2 = tst1 + 1.0_wp/tst1 
            if (tst2  >  tst1) go to 795 
            do 792 j = i, en 
               h(j,na) = h(j,na)/t 
               h(j,en) = h(j,en)/t 
  792       continue 
!                                                                       
  795    continue 
!     .......... end complex vector ..........                          
  800 continue 
!     .......... end back substitution.                                 
!                vectors of isolated roots ..........                   
      do 840 i = 1, n 
         if (i  >=  low .and. i  <=  igh) go to 840 
!                                                                       
         do 820 j = i, n 
  820    z(i,j) = h(i,j) 
!                                                                       
  840 continue 
!     .......... multiply by transformation matrix to give              
!                vectors of original full matrix.                       
!                for j=n step -1 until low do -- ..........             
      do 880 jj = low, n 
         j = n + low - jj 
         m = min0(j,igh) 
!                                                                       
         do 880 i = low, igh 
            zz = 0.0_wp 
!                                                                       
            do 860 k = low, m 
  860       zz = zz + z(i,k) * h(k,j) 
!                                                                       
            z(i,j) = zz 
  880 continue 
!                                                                       
      go to 1001 
!     .......... set error -- all eigenvalues have not                  
!                converged after 30*n iterations ..........             
 1000 ierr = en 
 1001 return 
      end subroutine hqr2

!=============================================================================80

      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr) 
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)   
!                                                                       
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr 
      real(wp) h(nm,n),wr(n),wi(n) 
      real(wp) p,q,r,s,t,w,x,y,zz,norm,tst1,tst2 
      logical notlas 
!                                                                       
!     this subroutine is a translation of the algol procedure hqr,      
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.    
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).   
!                                                                       
!     this subroutine finds the eigenvalues of a real                   
!     upper hessenberg matrix by the qr method.                         
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        low and igh are integers determined by the balancing           
!          subroutine  balanc.  if  balanc  has not been used,          
!          set low=1, igh=n.                                            
!                                                                       
!        h contains the upper hessenberg matrix.  information about     
!          the transformations used in the reduction to hessenberg      
!          form by  elmhes  or  orthes, if performed, is stored         
!          in the remaining triangle under the hessenberg matrix.       
!                                                                       
!     on output                                                         
!                                                                       
!        h has been destroyed.  therefore, it must be saved             
!          before calling  hqr  if subsequent calculation and           
!          back transformation of eigenvectors is to be performed.      
!                                                                       
!        wr and wi contain the real and imaginary parts,                
!          respectively, of the eigenvalues.  the eigenvalues           
!          are unordered except that complex conjugate pairs            
!          of values appear consecutively with the eigenvalue           
!          having the positive imaginary part first.  if an             
!          error exit is made, the eigenvalues should be correct        
!          for indices ierr+1,...,n.                                    
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the limit of 30*n iterations is exhausted      
!                     while the j-th eigenvalue is being sought.        
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated september 1989.                                
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      norm = 0.0_wp 
      k = 1 
!     .......... store roots isolated by balanc                         
!                and compute matrix norm ..........                     
      do 50 i = 1, n 
!                                                                       
         do 40 j = k, n 
   40    norm = norm + abs(h(i,j)) 
!                                                                       
         k = i 
         if (i  >=  low .and. i  <=  igh) go to 50 
         wr(i) = h(i,i) 
         wi(i) = 0.0_wp 
   50 continue 
!                                                                       
      en = igh 
      t = 0.0_wp 
      itn = 30*n 
!     .......... search for next eigenvalues ..........                 
   60 if (en  <  low) go to 1001 
      its = 0 
      na = en - 1 
      enm2 = na - 1 
!     .......... look for single small sub-diagonal element             
!                for l=en step -1 until low do -- ..........            
   70 do 80 ll = low, en 
         l = en + low - ll 
         if (l  ==  low) go to 100 
         s = abs(h(l-1,l-1)) + abs(h(l,l)) 
         if (s  ==  0.0_wp) s = norm 
         tst1 = s 
         tst2 = tst1 + abs(h(l,l-1)) 
         if (tst2  ==  tst1) go to 100 
   80 continue 
!     .......... form shift ..........                                  
  100 x = h(en,en) 
      if (l  ==  en) go to 270 
      y = h(na,na) 
      w = h(en,na) * h(na,en) 
      if (l  ==  na) go to 280 
      if (itn  ==  0) go to 1000 
      if (its  /=  10 .and. its  /=  20) go to 130 
!     .......... form exceptional shift ..........                      
      t = t + x 
!                                                                       
      do 120 i = low, en 
  120 h(i,i) = h(i,i) - x 
!                                                                       
      s = abs(h(en,na)) + abs(h(na,enm2)) 
      x = 0.75_wp * s 
      y = x 
      w = -0.4375_wp * s * s 
  130 its = its + 1 
      itn = itn - 1 
!     .......... look for two consecutive small                         
!                sub-diagonal elements.                                 
!                for m=en-2 step -1 until l do -- ..........            
      do 140 mm = l, enm2 
         m = enm2 + l - mm 
         zz = h(m,m) 
         r = x - zz 
         s = y - zz 
         p = (r * s - w) / h(m+1,m) + h(m,m+1) 
         q = h(m+1,m+1) - zz - r - s 
         r = h(m+2,m+1) 
         s = abs(p) + abs(q) + abs(r) 
         p = p / s 
         q = q / s 
         r = r / s 
         if (m  ==  l) go to 150 
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1))) 
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r)) 
         if (tst2  ==  tst1) go to 150 
  140 continue 
!                                                                       
  150 mp2 = m + 2 
!                                                                       
      do 160 i = mp2, en 
         h(i,i-2) = 0.0_wp 
         if (i  ==  mp2) go to 160 
         h(i,i-3) = 0.0_wp 
  160 continue 
!     .......... double qr step involving rows l to en and              
!                columns m to en ..........                             
      do 260 k = m, na 
         notlas = k  /=  na 
         if (k  ==  m) go to 170 
         p = h(k,k-1) 
         q = h(k+1,k-1) 
         r = 0.0_wp 
         if (notlas) r = h(k+2,k-1) 
         x = abs(p) + abs(q) + abs(r) 
         if (x  ==  0.0_wp) go to 260 
         p = p / x 
         q = q / x 
         r = r / x 
  170    s = sign(sqrt(p*p+q*q+r*r),p) 
         if (k  ==  m) go to 180 
         h(k,k-1) = -s * x 
         go to 190 
  180    if (l  /=  m) h(k,k-1) = -h(k,k-1) 
  190    p = p + s 
         x = p / s 
         y = q / s 
         zz = r / s 
         q = q / p 
         r = r / p 
         if (notlas) go to 225 
!     .......... row modification ..........                            
         do 200 j = k, EN 
            p = h(k,j) + q * h(k+1,j) 
            h(k,j) = h(k,j) - p * x 
            h(k+1,j) = h(k+1,j) - p * y 
  200    continue 
!                                                                       
         j = min0(en,k+3) 
!     .......... column modification ..........                         
         do 210 i = L, j 
            p = x * h(i,k) + y * h(i,k+1) 
            h(i,k) = h(i,k) - p 
            h(i,k+1) = h(i,k+1) - p * q 
  210    continue 
         go to 255 
  225    continue 
!     .......... row modification ..........                            
         do 230 j = k, EN 
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j) 
            h(k,j) = h(k,j) - p * x 
            h(k+1,j) = h(k+1,j) - p * y 
            h(k+2,j) = h(k+2,j) - p * zz 
  230    continue 
!                                                                       
         j = min0(en,k+3) 
!     .......... column modification ..........                         
         do 240 i = L, j 
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2) 
            h(i,k) = h(i,k) - p 
            h(i,k+1) = h(i,k+1) - p * q 
            h(i,k+2) = h(i,k+2) - p * r 
  240    continue 
  255    continue 
!                                                                       
  260 continue 
!                                                                       
      go to 70 
!     .......... one root found ..........                              
  270 wr(en) = x + t 
      wi(en) = 0.0_wp 
      en = na 
      go to 60 
!     .......... two roots found ..........                             
  280 p = (y - x) / 2.0_wp 
      q = p * p + w 
      zz = sqrt(abs(q)) 
      x = x + t 
      if (q  <  0.0_wp) go to 320 
!     .......... real pair ..........                                   
      zz = p + sign(zz,p) 
      wr(na) = x + zz 
      wr(en) = wr(na) 
      if (zz  /=  0.0_wp) wr(en) = x - w / zz 
      wi(na) = 0.0_wp 
      wi(en) = 0.0_wp 
      go to 330 
!     .......... complex pair ..........                                
  320 wr(na) = x + p 
      wr(en) = x + p 
      wi(na) = zz 
      wi(en) = -zz 
  330 en = enm2 
      go to 60 
!     .......... set error -- all eigenvalues have not                  
!                converged after 30*n iterations ..........             
 1000 ierr = en 
 1001 return 
      end subroutine hqr

!=============================================================================80

      subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1) 
!                                                                       
      integer i,j,k,l,m,n,ii,mml,tag,ierr 
      real(wp) d(n),e(n),e2(n),w(n),rv1(n) 
      real(wp) b,c,f,g,p,r,s,tst1,tst2
      integer ind(n) 
!                                                                       
!     this subroutine is a variant of  imtql1  which is a translation of
!     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
!     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.   
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).   
!                                                                       
!     this subroutine finds the eigenvalues of a symmetric tridiagonal  
!     matrix by the implicit ql method and associates with them         
!     their corresponding submatrix indices.                            
!                                                                       
!     on input                                                          
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        d contains the diagonal elements of the input matrix.          
!                                                                       
!        e contains the subdiagonal elements of the input matrix        
!          in its last n-1 positions.  e(1) is arbitrary.               
!                                                                       
!        e2 contains the squares of the corresponding elements of e.    
!          e2(1) is arbitrary.                                          
!                                                                       
!     on output                                                         
!                                                                       
!        d and e are unaltered.                                         
!                                                                       
!        elements of e2, corresponding to elements of e regarded        
!          as negligible, have been replaced by zero causing the        
!          matrix to split into a direct sum of submatrices.            
!          e2(1) is also set to zero.                                   
!                                                                       
!        w contains the eigenvalues in ascending order.  if an          
!          error exit is made, the eigenvalues are correct and          
!          ordered for indices 1,2,...ierr-1, but may not be            
!          the smallest eigenvalues.                                    
!                                                                       
!        ind contains the submatrix indices associated with the         
!          corresponding eigenvalues in w -- 1 for eigenvalues          
!          belonging to the first submatrix from the top,               
!          2 for those belonging to the second submatrix, etc..         
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the j-th eigenvalue has not been               
!                     determined after 30 iterations.                   
!                                                                       
!        rv1 is a temporary storage array.                              
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      k = 0 
      tag = 0 
!                                                                       
      do 100 i = 1, n 
         w(i) = d(i) 
         if (i  /=  1) rv1(i-1) = e(i) 
  100 continue 
!                                                                       
      e2(1) = 0.0_wp 
      rv1(n) = 0.0_wp 
!                                                                       
      do 290 l = 1, n 
         j = 0 
!     .......... look for small sub-diagonal element ..........         
  105    do 110 m = l, n 
            if (m  ==  n) go to 120 
            tst1 = abs(w(m)) + abs(w(m+1)) 
            tst2 = tst1 + abs(rv1(m)) 
            if (tst2  ==  tst1) go to 120 
!     .......... guard against underflowed element of e2 ..........     
            if (e2(m+1)  ==  0.0_wp) go to 125 
  110    continue 
!                                                                       
  120    if (m  <=  k) go to 130 
         if (m  /=  n) e2(m+1) = 0.0_wp 
  125    k = m 
         tag = tag + 1 
  130    p = w(l) 
         if (m  ==  l) go to 215 
         if (j  ==  30) go to 1000 
         j = j + 1 
!     .......... form shift ..........                                  
         g = (w(l+1) - p) / (2.0_wp * rv1(l)) 
         r = pythag(g,1.0_wp) 
         g = w(m) - p + rv1(l) / (g + sign(r,g)) 
         s = 1.0_wp 
         c = 1.0_wp 
         p = 0.0_wp 
         mml = m - l 
!     .......... for i=m-1 step -1 until l do -- ..........             
         do 200 ii = 1, mml 
            i = m - ii 
            f = s * rv1(i) 
            b = c * rv1(i) 
            r = pythag(f,g) 
            rv1(i+1) = r 
            if (r  ==  0.0_wp) go to 210 
            s = f / r 
            c = g / r 
            g = w(i+1) - p 
            r = (w(i) - g) * s + 2.0_wp * c * b 
            p = s * r 
            w(i+1) = g + p 
            g = c * r - b 
  200    continue 
!                                                                       
         w(l) = w(l) - p 
         rv1(l) = g 
         rv1(m) = 0.0_wp 
         go to 105 
!     .......... recover from underflow ..........                      
  210    w(i+1) = w(i+1) - p 
         rv1(m) = 0.0_wp 
         go to 105 
!     .......... order eigenvalues ..........                           
  215    if (l  ==  1) go to 250 
!     .......... for i=l step -1 until 2 do -- ..........               
         do 230 ii = 2, l 
            i = l + 2 - ii 
            if (p  >=  w(i-1)) go to 270 
            w(i) = w(i-1) 
            ind(i) = ind(i-1) 
  230    continue 
!                                                                       
  250    i = 1 
  270    w(i) = p 
         ind(i) = tag 
  290 continue 
!                                                                       
      go to 1001 
!     .......... set error -- no convergence to an                      
!                eigenvalue after 30 iterations ..........              
 1000 ierr = l 
 1001 return 
      end subroutine imtqlv

!=============================================================================80

      subroutine minfit(nm,m,n,a,w,ip,b,ierr,rv1) 
!                                                                       
      integer i,j,k,l,m,n,ii,ip,i1,kk,k1,ll,l1,m1,nm,its,ierr 
      real(wp) a(nm,n),w(n),b(nm,ip),rv1(n) 
      real(wp) c,f,g,h,s,x,y,z,tst1,tst2,scale
!                                                                       
!     this subroutine is a translation of the algol procedure minfit,   
!     num. math. 14, 403-420(1970) by golub and reinsch.                
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).   
!                                                                       
!     this subroutine determines, towards the solution of the linear    
!                                                        t              
!     system ax=b, the singular value decomposition a=usv  of a real    
!                                         t                             
!     m by n rectangular matrix, forming u b rather than u.  householder
!     bidiagonalization and a variant of the qr algorithm are used.     
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.  note that nm must be at least          
!          as large as the maximum of m and n.                          
!                                                                       
!        m is the number of rows of a and b.                            
!                                                                       
!        n is the number of columns of a and the order of v.            
!                                                                       
!        a contains the rectangular coefficient matrix of the system.   
!                                                                       
!        ip is the number of columns of b.  ip can be zero.             
!                                                                       
!        b contains the constant column matrix of the system            
!          if ip is not zero.  otherwise b is not referenced.           
!                                                                       
!     on output                                                         
!                                                                       
!        a has been overwritten by the matrix v (orthogonal) of the     
!          decomposition in its first n rows and columns.  if an        
!          error exit is made, the columns of v corresponding to        
!          indices of correct singular values should be correct.        
!                                                                       
!        w contains the n (non-negative) singular values of a (the      
!          diagonal elements of s).  they are unordered.  if an         
!          error exit is made, the singular values should be correct    
!          for indices ierr+1,ierr+2,...,n.                             
!                                                                       
!                                   t                                   
!        b has been overwritten by u b.  if an error exit is made,      
!                       t                                               
!          the rows of u b corresponding to indices of correct          
!          singular values should be correct.                           
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          k          if the k-th singular value has not been           
!                     determined after 30 iterations.                   
!                                                                       
!        rv1 is a temporary storage array.                              
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
!     .......... householder reduction to bidiagonal form ..........    
      g = 0.0_wp 
      scale = 0.0_wp 
      x = 0.0_wp 
!                                                                       
      do 300 i = 1, n 
         l = i + 1 
         rv1(i) = scale * g 
         g = 0.0_wp 
         s = 0.0_wp 
         scale = 0.0_wp 
         if (i  >  m) go to 210 
!                                                                       
         do 120 k = i, m 
  120    scale = scale + abs(a(k,i)) 
!                                                                       
         if (scale  ==  0.0_wp) go to 210 
!                                                                       
         do 130 k = i, m 
            a(k,i) = a(k,i) / scale 
            s = s + a(k,i)**2 
  130    continue 
!                                                                       
         f = a(i,i) 
         g = -sign(sqrt(s),f) 
         h = f * g - s 
         a(i,i) = f - g 
         if (i  ==  n) go to 160 
!                                                                       
         do 150 j = l, n 
            s = 0.0_wp 
!                                                                       
            do 140 k = i, m 
  140       s = s + a(k,i) * a(k,j) 
!                                                                       
            f = s / h 
!                                                                       
            do 150 k = i, m 
               a(k,j) = a(k,j) + f * a(k,i) 
  150    continue 
!                                                                       
  160    if (ip  ==  0) go to 190 
!                                                                       
         do 180 j = 1, ip 
            s = 0.0_wp 
!                                                                       
            do 170 k = i, m 
  170       s = s + a(k,i) * b(k,j) 
!                                                                       
            f = s / h 
!                                                                       
            do 180 k = i, m 
               b(k,j) = b(k,j) + f * a(k,i) 
  180    continue 
!                                                                       
  190    do 200 k = i, m 
  200    a(k,i) = scale * a(k,i) 
!                                                                       
  210    w(i) = scale * g 
         g = 0.0_wp 
         s = 0.0_wp 
         scale = 0.0_wp 
         if (i  >  m .or. i  ==  n) go to 290 
!                                                                       
         do 220 k = l, n 
  220    scale = scale + abs(a(i,k)) 
!                                                                       
         if (scale  ==  0.0_wp) go to 290 
!                                                                       
         do 230 k = l, n 
            a(i,k) = a(i,k) / scale 
            s = s + a(i,k)**2 
  230    continue 
!                                                                       
         f = a(i,l) 
         g = -sign(sqrt(s),f) 
         h = f * g - s 
         a(i,l) = f - g 
!                                                                       
         do 240 k = l, n 
  240    rv1(k) = a(i,k) / h 
!                                                                       
         if (i  ==  m) go to 270 
!                                                                       
         do 260 j = l, m 
            s = 0.0_wp 
!                                                                       
            do 250 k = l, n 
  250       s = s + a(j,k) * a(i,k) 
!                                                                       
            do 260 k = l, n 
               a(j,k) = a(j,k) + s * rv1(k) 
  260    continue 
!                                                                       
  270    do 280 k = l, n 
  280    a(i,k) = scale * a(i,k) 
!                                                                       
  290    x = max(x,abs(w(i))+abs(rv1(i))) 
  300 continue 
!     .......... accumulation of right-hand transformations.            
!                for i=n step -1 until 1 do -- ..........               
      do 400 ii = 1, n 
         i = n + 1 - ii 
         if (i  ==  n) go to 390 
         if (g  ==  0.0_wp) go to 360 
!                                                                       
         do 320 j = l, n 
!     .......... double division avoids possible underflow ..........   
  320    a(j,i) = (a(i,j) / a(i,l)) / g 
!                                                                       
         do 350 j = l, n 
            s = 0.0_wp 
!                                                                       
            do 340 k = l, n 
  340       s = s + a(i,k) * a(k,j) 
!                                                                       
            do 350 k = l, n 
               a(k,j) = a(k,j) + s * a(k,i) 
  350    continue 
!                                                                       
  360    do 380 j = l, n 
            a(i,j) = 0.0_wp 
            a(j,i) = 0.0_wp 
  380    continue 
!                                                                       
  390    a(i,i) = 1.0_wp 
         g = rv1(i) 
         l = i 
  400 continue 
!                                                                       
      if (m  >=  n .or. ip  ==  0) go to 510 
      m1 = m + 1 
!                                                                       
      do 500 i = m1, n 
!                                                                       
         do 500 j = 1, ip 
            b(i,j) = 0.0_wp 
  500 continue 
!     .......... diagonalization of the bidiagonal form ..........      
  510 tst1 = x 
!     .......... for k=n step -1 until 1 do -- ..........               
      do 700 kk = 1, n 
         k1 = n - kk 
         k = k1 + 1 
         its = 0 
!     .......... test for splitting.                                    
!                for l=k step -1 until 1 do -- ..........               
  520    do 530 ll = 1, k 
            l1 = k - ll 
            l = l1 + 1 
            tst2 = tst1 + abs(rv1(l)) 
            if (tst2  ==  tst1) go to 565 
!     .......... rv1(1) is always zero, so there is no exit             
!                through the bottom of the loop ..........              
            tst2 = tst1 + abs(w(l1)) 
            if (tst2  ==  tst1) go to 540 
  530    continue 
!     .......... cancellation of rv1(l) if l greater than 1 ..........  
  540    c = 0.0_wp 
         s = 1.0_wp 
!                                                                       
         do 560 i = l, k 
            f = s * rv1(i) 
            rv1(i) = c * rv1(i) 
            tst2 = tst1 + abs(f) 
            if (tst2  ==  tst1) go to 565 
            g = w(i) 
            h = pythag(f,g) 
            w(i) = h 
            c = g / h 
            s = -f / h 
            if (ip  ==  0) go to 560 
!                                                                       
            do 550 j = 1, ip 
               y = b(l1,j) 
               z = b(i,j) 
               b(l1,j) = y * c + z * s 
               b(i,j) = -y * s + z * c 
  550       continue 
!                                                                       
  560    continue 
!     .......... test for convergence ..........                        
  565    z = w(k) 
         if (l  ==  k) go to 650 
!     .......... shift from bottom 2 by 2 minor ..........              
         if (its  ==  30) go to 1000 
         its = its + 1 
         x = w(l) 
         y = w(k1) 
         g = rv1(k1) 
         h = rv1(k) 
         f = 0.5_wp* (((g + z) / h) * ((g - z) / y) + y / h - h / y) 
         g = pythag(f,1.0_wp) 
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h) 
!     .......... next qr transformation ..........                      
         c = 1.0_wp 
         s = 1.0_wp 
!                                                                       
         do 600 i1 = l, k1 
            i = i1 + 1 
            g = rv1(i) 
            y = w(i) 
            h = s * g 
            g = c * g 
            z = pythag(f,h) 
            rv1(i1) = z 
            c = f / z 
            s = h / z 
            f = x * c + g * s 
            g = -x * s + g * c 
            h = y * s 
            y = y * c 
!                                                                       
            do 570 j = 1, n 
               x = a(j,i1) 
               z = a(j,i) 
               a(j,i1) = x * c + z * s 
               a(j,i) = -x * s + z * c 
  570       continue 
!                                                                       
            z = pythag(f,h) 
            w(i1) = z 
!     .......... rotation can be arbitrary if z is zero ..........      
            if (z  ==  0.0_wp) go to 580 
            c = f / z 
            s = h / z 
  580       f = c * g + s * y 
            x = -s * g + c * y 
            if (ip  ==  0) go to 600 
!                                                                       
            do 590 j = 1, ip 
               y = b(i1,j) 
               z = b(i,j) 
               b(i1,j) = y * c + z * s 
               b(i,j) = -y * s + z * c 
  590       continue 
!                                                                       
  600    continue 
!                                                                       
         rv1(l) = 0.0_wp 
         rv1(k) = f 
         w(k) = x 
         go to 520 
!     .......... convergence ..........                                 
  650    if (z  >=  0.0_wp) go to 700 
!     .......... w(k) is made non-negative ..........                   
         w(k) = -z 
!                                                                       
         do 690 j = 1, n 
  690    a(j,k) = -a(j,k) 
!                                                                       
  700 continue 
!                                                                       
      go to 1001 
!     .......... set error -- no convergence to a                       
!                singular value after 30 iterations ..........          
 1000 ierr = k 
 1001 return 
      end subroutine minfit

!=============================================================================80

      subroutine tqlrat(n,d,e2,ierr) 
!                                                                       
      integer i,j,l,m,n,ii,l1,mml,ierr 
      real(wp) d(n),e2(n) 
      real(wp) b,c,f,g,h,p,r,s,t
!                                                                       
!     this subroutine is a translation of the algol procedure tqlrat,   
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.                
!                                                                       
!     this subroutine finds the eigenvalues of a symmetric              
!     tridiagonal matrix by the rational ql method.                     
!                                                                       
!     on input                                                          
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        d contains the diagonal elements of the input matrix.          
!                                                                       
!        e2 contains the squares of the subdiagonal elements of the     
!          input matrix in its last n-1 positions.  e2(1) is arbitrary. 
!                                                                       
!      on output                                                        
!                                                                       
!        d contains the eigenvalues in ascending order.  if an          
!          error exit is made, the eigenvalues are correct and          
!          ordered for indices 1,2,...ierr-1, but may not be            
!          the smallest eigenvalues.                                    
!                                                                       
!        e2 has been destroyed.                                         
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the j-th eigenvalue has not been               
!                     determined after 30 iterations.                   
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      if (n  ==  1) go to 1001 
!                                                                       
      do 100 i = 2, n 
  100 e2(i-1) = e2(i) 
!                                                                       
      f = 0.0_wp 
      t = 0.0_wp 
      e2(n) = 0.0_wp 
!                                                                       
      do 290 l = 1, n 
         j = 0 
         h = abs(d(l)) + sqrt(e2(l)) 
         if (t  >  h) go to 105 
         t = h 
         b = epslon(t) 
         c = b * b 
!     .......... look for small squared sub-diagonal element .......... 
  105    do 110 m = l, n 
            if (e2(m)  <=  c) go to 120 
!     .......... e2(n) is always zero, so there is no exit              
!                through the bottom of the loop ..........              
  110    continue 
!                                                                       
  120    if (m  ==  l) go to 210 
  130    if (j  ==  50) go to 1000 
         j = j + 1 
!     .......... form shift ..........                                  
         l1 = l + 1 
         s = sqrt(e2(l)) 
         g = d(l) 
         p = (d(l1) - g) / (2.0_wp * s) 
         r = pythag(p,1.0_wp) 
         d(l) = s / (p + sign(r,p)) 
         h = g - d(l) 
!                                                                       
         do 140 i = l1, n 
  140    d(i) = d(i) - h 
!                                                                       
         f = f + h 
!     .......... rational ql transformation ..........                  
         g = d(m) 
         if (g  ==  0.0_wp) g = b 
         h = g 
         s = 0.0_wp 
         mml = m - l 
!     .......... for i=m-1 step -1 until l do -- ..........             
         do 200 ii = 1, mml 
            i = m - ii 
            p = g * h 
            r = p + e2(i) 
            e2(i+1) = s * r 
            s = e2(i) / r 
            d(i+1) = h + s * (h + d(i)) 
            g = d(i) - e2(i) / g 
            if (g  ==  0.0_wp) g = b 
            h = g * p / r 
  200    continue 
!                                                                       
         e2(l) = s * g 
         d(l) = h 
!     .......... guard against underflow in convergence test .......... 
         if (h  ==  0.0_wp) go to 210 
         if (abs(e2(l))  <=  abs(c/h)) go to 210 
         e2(l) = h * e2(l) 
         if (e2(l)  /=  0.0_wp) go to 130 
  210    p = d(l) + f 
!     .......... order eigenvalues ..........                           
         if (l  ==  1) go to 250 
!     .......... for i=l step -1 until 2 do -- ..........               
         do 230 ii = 2, l 
            i = l + 2 - ii 
            if (p  >=  d(i-1)) go to 270 
            d(i) = d(i-1) 
  230    continue 
!                                                                       
  250    i = 1 
  270    d(i) = p 
  290 continue 
!                                                                       
      go to 1001 
!     .......... set error -- no convergence to an                      
!                eigenvalue after 30 iterations ..........              
 1000 ierr = l 
 1001 return 
      end subroutine tqlrat

!=============================================================================80

      real(wp) function pythag(a,b) 
      real(wp) a,b 
!                                                                       
!     finds sqrt(a**2+b**2) without overflow or destructive underflow  
!                                                                       
      real(wp) p,r,s,t,u 
      p = max(abs(a),abs(b)) 
      if (p  ==  0.0_wp) go to 20 
      r = (min(abs(a),abs(b))/p)**2 
   10 continue 
         t = 4.0_wp + r 
         if (t  ==  4.0_wp) go to 20 
         s = r/t 
         u = 1.0_wp + 2.0_wp*s 
         p = u*p 
         r = (s/u)**2 * r 
      go to 10 
   20 pythag = p 
      return 
      end function pythag

!=============================================================================80

      subroutine qsortd(x,ind,n)

!       Code converted using TO_F90 by Alan Miller
!       Date: 2002-12-18  Time: 11:55:47

      real(wp), INTENT(IN)  :: x(:)
      INTEGER, INTENT(OUT)   :: ind(:)
      INTEGER, INTENT(IN)    :: n

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
      REAL(wp) :: t

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

      subroutine qzhes(nm,n,a,b,matz,z) 
!                                                                       
!     integer, parameter               :: dp = 8    ! double precision
      integer i,j,k,l,n,lb,l1,nm,nk1,nm1,nm2 
      real(wp) a(nm,n),b(nm,n),z(nm,n) 
      real(wp) r,s,t,u1,u2,v1,v2,rho 
      logical matz 
!                                                                       
!     this subroutine is the first step of the qz algorithm             
!     for solving generalized matrix eigenvalue problems,               
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.      
!                                                                       
!     this subroutine accepts a pair of real general matrices and       
!     reduces one of them to upper hessenberg form and the other        
!     to upper triangular form using orthogonal transformations.        
!     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.   
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrices.                                
!                                                                       
!        a contains a real general matrix.                              
!                                                                       
!        b contains a real general matrix.                              
!                                                                       
!        matz should be set to .true. if the right hand transformations 
!          are to be accumulated for later use in computing             
!          eigenvectors, and to .false. otherwise.                      
!                                                                       
!     on output                                                         
!                                                                       
!        a has been reduced to upper hessenberg form.  the elements     
!          below the first subdiagonal have been set to zero.           
!                                                                       
!        b has been reduced to upper triangular form.  the elements     
!          below the main diagonal have been set to zero.               
!                                                                       
!        z contains the product of the right hand transformations if    
!          matz has been set to .true.  otherwise, z is not referenced. 
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
!     .......... initialize z ..........                                
      if (.not. matz) go to 10 
!                                                                       
      do 3 j = 1, n 
!                                                                       
         do 2 i = 1, n 
            z(i,j) = 0.0_wp 
    2    continue 
!                                                                       
         z(j,j) = 1.0_wp 
    3 continue 
!     .......... reduce b to upper triangular form ..........           
   10 if (n  <=  1) go to 170 
      nm1 = n - 1 
!                                                                       
      do 100 l = 1, nm1 
         l1 = l + 1 
         s = 0.0_wp 
!                                                                       
         do 20 i = l1, n 
            s = s + abs(b(i,l)) 
   20    continue 
!                                                                       
         if (s  ==  0.0_wp) go to 100 
         s = s + abs(b(l,l)) 
         r = 0.0_wp 
!                                                                       
         do 25 i = l, n 
            b(i,l) = b(i,l) / s 
            r = r + b(i,l)**2 
   25    continue 
!                                                                       
         r = sign(sqrt(r),b(l,l)) 
         b(l,l) = b(l,l) + r 
         rho = r * b(l,l) 
!                                                                       
         do 50 j = l1, n 
            t = 0.0_wp 
!                                                                       
            do 30 i = l, n 
               t = t + b(i,l) * b(i,j) 
   30       continue 
!                                                                       
            t = -t / rho 
!                                                                       
            do 40 i = l, n 
               b(i,j) = b(i,j) + t * b(i,l) 
   40       continue 
!                                                                       
   50    continue 
!                                                                       
         do 80 j = 1, n 
            t = 0.0_wp 
!                                                                       
            do 60 i = l, n 
               t = t + b(i,l) * a(i,j) 
   60       continue 
!                                                                       
            t = -t / rho 
!                                                                       
            do 70 i = l, n 
               a(i,j) = a(i,j) + t * b(i,l) 
   70       continue 
!                                                                       
   80    continue 
!                                                                       
         b(l,l) = -s * r 
!                                                                       
         do 90 i = l1, n 
            b(i,l) = 0.0_wp 
   90    continue 
!                                                                       
  100 continue 
!     .......... reduce a to upper hessenberg form, while               
!                keeping b triangular ..........                        
      if (n  ==  2) go to 170 
      nm2 = n - 2 
!                                                                       
      do 160 k = 1, nm2 
         nk1 = nm1 - k 
!     .......... for l=n-1 step -1 until k+1 do -- ..........           
         do 150 lb = 1, nk1 
            l = n - lb 
            l1 = l + 1 
!     .......... zero a(l+1,k) ..........                               
            s = abs(a(l,k)) + abs(a(l1,k)) 
            if (s  ==  0.0_wp) go to 150 
            u1 = a(l,k) / s 
            u2 = a(l1,k) / s 
            r = sign(sqrt(u1*u1+u2*u2),u1) 
            v1 =  -(u1 + r) / r 
            v2 = -u2 / r 
            u2 = v2 / v1 
!                                                                       
            do 110 j = k, n 
               t = a(l,j) + u2 * a(l1,j) 
               a(l,j) = a(l,j) + t * v1 
               a(l1,j) = a(l1,j) + t * v2 
  110       continue 
!                                                                       
            a(l1,k) = 0.0_wp 
!                                                                       
            do 120 j = l, n 
               t = b(l,j) + u2 * b(l1,j) 
               b(l,j) = b(l,j) + t * v1 
               b(l1,j) = b(l1,j) + t * v2 
  120       continue 
!     .......... zero b(l+1,l) ..........                               
            s = abs(b(l1,l1)) + abs(b(l1,l)) 
            if (s  ==  0.0_wp) go to 150 
            u1 = b(l1,l1) / s 
            u2 = b(l1,l) / s 
            r = sign(sqrt(u1*u1+u2*u2),u1) 
            v1 =  -(u1 + r) / r 
            v2 = -u2 / r 
            u2 = v2 / v1 
!                                                                       
            do 130 i = 1, l1 
               t = b(i,l1) + u2 * b(i,l) 
               b(i,l1) = b(i,l1) + t * v1 
               b(i,l) = b(i,l) + t * v2 
  130       continue 
!                                                                       
            b(l1,l) = 0.0_wp 
!                                                                       
            do 140 i = 1, n 
               t = a(i,l1) + u2 * a(i,l) 
               a(i,l1) = a(i,l1) + t * v1 
               a(i,l) = a(i,l) + t * v2 
  140       continue 
!                                                                       
            if (.not. matz) go to 150 
!                                                                       
            do 145 i = 1, n 
               t = z(i,l1) + u2 * z(i,l) 
               z(i,l1) = z(i,l1) + t * v1 
               z(i,l) = z(i,l) + t * v2 
  145       continue 
!                                                                       
  150    continue 
!                                                                       
  160 continue 
!                                                                       
  170 return 
      end subroutine qzhes

!=============================================================================80

      subroutine qzit(nm,n,a,b,eps1,matz,z,ierr) 
!                                                                       
!     integer, parameter               :: dp = 8    ! double precision
      integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,itn,its,km1,lm1,    &
     &        enm2,ierr,lor1,enorn                                      
      real(wp) a(nm,n),b(nm,n),z(nm,n) 
      real(wp) r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11,  &
     &       a12,a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34,       &
     &       b44,epsa,epsb,eps1,anorm,bnorm
      logical matz,notlas 
!                                                                       
!     this subroutine is the second step of the qz algorithm            
!     for solving generalized matrix eigenvalue problems,               
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,      
!     as modified in technical note nasa tn d-7305(1973) by ward.       
!                                                                       
!     this subroutine accepts a pair of real matrices, one of them      
!     in upper hessenberg form and the other in upper triangular form.  
!     it reduces the hessenberg matrix to quasi-triangular form using   
!     orthogonal transformations while maintaining the triangular form  
!     of the other matrix.  it is usually preceded by  qzhes  and       
!     followed by  qzval  and, possibly,  qzvec.                        
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrices.                                
!                                                                       
!        a contains a real upper hessenberg matrix.                     
!                                                                       
!        b contains a real upper triangular matrix.                     
!                                                                       
!        eps1 is a tolerance used to determine negligible elements.     
!          eps1 = 0.0 (or negative) may be input, in which case an      
!          element will be neglected only if it is less than roundoff   
!          error times the norm of its matrix.  if the input eps1 is    
!          positive, then an element will be considered negligible      
!          if it is less than eps1 times the norm of its matrix.  a     
!          positive value of eps1 may result in faster execution,       
!          but less accurate results.                                   
!                                                                       
!        matz should be set to .true. if the right hand transformations 
!          are to be accumulated for later use in computing             
!          eigenvectors, and to .false. otherwise.                      
!                                                                       
!        z contains, if matz has been set to .true., the                
!          transformation matrix produced in the reduction              
!          by  qzhes, if performed, or else the identity matrix.        
!          if matz has been set to .false., z is not referenced.        
!                                                                       
!     on output                                                         
!                                                                       
!        a has been reduced to quasi-triangular form.  the elements     
!          below the first subdiagonal are still zero and no two        
!          consecutive subdiagonal elements are nonzero.                
!                                                                       
!        b is still in upper triangular form, although its elements     
!          have been altered.  the location b(n,1) is used to store     
!          eps1 times the norm of b for later use by  qzval  and  qzvec.
!                                                                       
!        z contains the product of the right hand transformations       
!          (for both steps) if matz has been set to .true..             
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the limit of 30*n iterations is exhausted      
!                     while the j-th eigenvalue is being sought.        
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
!     .......... compute epsa,epsb ..........                           
      anorm = 0.0_wp 
      bnorm = 0.0_wp 
!                                                                       
      do 30 i = 1, n 
         ani = 0.0_wp 
         if (i  /=  1) ani = abs(a(i,i-1)) 
         bni = 0.0_wp 
!                                                                       
         do 20 j = i, n 
            ani = ani + abs(a(i,j)) 
            bni = bni + abs(b(i,j)) 
   20    continue 
!                                                                       
         if (ani  >  anorm) anorm = ani 
         if (bni  >  bnorm) bnorm = bni 
   30 continue 
!                                                                       
      if (anorm  ==  0.0_wp) anorm = 1.0_wp 
      if (bnorm  ==  0.0_wp) bnorm = 1.0_wp 
      ep = eps1 
      if (ep  >  0.0_wp) go to 50 
!     .......... use roundoff level if eps1 is zero ..........          
      ep = epslon(1.0_wp) 
   50 epsa = ep * anorm 
      epsb = ep * bnorm 
!     .......... reduce a to quasi-triangular form, while               
!                keeping b triangular ..........                        
      lor1 = 1 
      enorn = n 
      en = n 
!     itn = 30*n   original iteration count
      itn = 50*n 
!     .......... begin qz step ..........                               
   60 if (en  <=  2) go to 1001 
      if (.not. matz) enorn = en 
      its = 0 
      na = en - 1 
      enm2 = na - 1 
   70 ish = 2 
!     .......... check for convergence or reducibility.                 
!                for l=en step -1 until 1 do -- ..........              
      do 80 ll = 1, en 
         lm1 = en - ll 
         l = lm1 + 1 
         if (l  ==  1) go to 95 
         if (abs(a(l,lm1))  <=  epsa) go to 90 
   80 continue 
!                                                                       
   90 a(l,lm1) = 0.0_wp 
      if (l  <  na) go to 95 
!     .......... 1-by-1 or 2-by-2 block isolated ..........             
      en = lm1 
      go to 60 
!     .......... check for small top of b ..........                    
   95 ld = l 
  100 l1 = l + 1 
      b11 = b(l,l) 
      if (abs(b11)  >  epsb) go to 120 
      b(l,l) = 0.0_wp 
      s = abs(a(l,l)) + abs(a(l1,l)) 
      u1 = a(l,l) / s 
      u2 = a(l1,l) / s 
      r = sign(sqrt(u1*u1+u2*u2),u1) 
      v1 = -(u1 + r) / r 
      v2 = -u2 / r 
      u2 = v2 / v1 
!                                                                       
      do 110 j = l, enorn 
         t = a(l,j) + u2 * a(l1,j) 
         a(l,j) = a(l,j) + t * v1 
         a(l1,j) = a(l1,j) + t * v2 
         t = b(l,j) + u2 * b(l1,j) 
         b(l,j) = b(l,j) + t * v1 
         b(l1,j) = b(l1,j) + t * v2 
  110 continue 
!                                                                       
      if (l  /=  1) a(l,lm1) = -a(l,lm1) 
      lm1 = l 
      l = l1 
      go to 90 
  120 a11 = a(l,l) / b11 
      a21 = a(l1,l) / b11 
      if (ish  ==  1) go to 140 
!     .......... iteration strategy ..........                          
      if (itn  ==  0) go to 1000 
      if (its  ==  10) go to 155 
!     .......... determine type of shift ..........                     
      b22 = b(l1,l1) 
      if (abs(b22)  <  epsb) b22 = epsb 
      b33 = b(na,na) 
      if (abs(b33)  <  epsb) b33 = epsb 
      b44 = b(en,en) 
      if (abs(b44)  <  epsb) b44 = epsb 
      a33 = a(na,na) / b33 
      a34 = a(na,en) / b44 
      a43 = a(en,na) / b33 
      a44 = a(en,en) / b44 
      b34 = b(na,en) / b44 
      t = 0.5_wp * (a43 * b34 - a33 - a44) 
      r = t * t + a34 * a43 - a33 * a44 
      if (r  <  0.0_wp) go to 150 
!     .......... determine single shift zeroth column of a ..........   
      ish = 1 
      r = sqrt(r) 
      sh = -t + r 
      s = -t - r 
      if (abs(s-a44)  <  abs(sh-a44)) sh = s 
!     .......... look for two consecutive small                         
!                sub-diagonal elements of a.                            
!                for l=en-2 step -1 until ld do -- ..........           
      do 130 ll = ld, enm2 
         l = enm2 + ld - ll 
         if (l  ==  ld) go to 140 
         lm1 = l - 1 
         l1 = l + 1 
         t = a(l,l) 
         if (abs(b(l,l))  >  epsb) t = t - sh * b(l,l) 
         if (abs(a(l,lm1))  <=  abs(t/a(l1,l)) * epsa) go to 100 
  130 continue 
!                                                                       
  140 a1 = a11 - sh 
      a2 = a21 
      if (l  /=  ld) a(l,lm1) = -a(l,lm1) 
      go to 160 
!     .......... determine double shift zeroth column of a ..........   
  150 a12 = a(l,l1) / b22 
      a22 = a(l1,l1) / b22 
      b12 = b(l,l1) / b22 
      a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11)    &
     &     / a21 + a12 - a11 * b12                                      
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11)          &
     &     + a43 * b34                                                  
      a3 = a(l1+1,l1) / b22 
      go to 160 
!     .......... ad hoc shift ..........                                
  155 a1 = 0.0_wp 
      a2 = 1.0_wp 
      a3 = 1.1605_wp 
  160 its = its + 1 
      itn = itn - 1 
      if (.not. matz) lor1 = ld 
!     .......... main loop ..........                                   
      do 260 k = l, na 
         notlas = k  /=  na .and. ish  ==  2 
         k1 = k + 1 
         k2 = k + 2 
         km1 = max0(k-1,l) 
         ll = min0(en,k1+ish) 
         if (notlas) go to 190 
!     .......... zero a(k+1,k-1) ..........                             
         if (k  ==  l) go to 170 
         a1 = a(k,km1) 
         a2 = a(k1,km1) 
  170    s = abs(a1) + abs(a2) 
         if (s  ==  0.0_wp) go to 70 
         u1 = a1 / s 
         u2 = a2 / s 
         r = sign(sqrt(u1*u1+u2*u2),u1) 
         v1 = -(u1 + r) / r 
         v2 = -u2 / r 
         u2 = v2 / v1 
!                                                                       
         do 180 j = km1, enorn 
            t = a(k,j) + u2 * a(k1,j) 
            a(k,j) = a(k,j) + t * v1 
            a(k1,j) = a(k1,j) + t * v2 
            t = b(k,j) + u2 * b(k1,j) 
            b(k,j) = b(k,j) + t * v1 
            b(k1,j) = b(k1,j) + t * v2 
  180    continue 
!                                                                       
         if (k  /=  l) a(k1,km1) = 0.0_wp 
         go to 240 
!     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........              
  190    if (k  ==  l) go to 200 
         a1 = a(k,km1) 
         a2 = a(k1,km1) 
         a3 = a(k2,km1) 
  200    s = abs(a1) + abs(a2) + abs(a3) 
         if (s  ==  0.0_wp) go to 260 
         u1 = a1 / s 
         u2 = a2 / s 
         u3 = a3 / s 
         r = sign(sqrt(u1*u1+u2*u2+u3*u3),u1) 
         v1 = -(u1 + r) / r 
         v2 = -u2 / r 
         v3 = -u3 / r 
         u2 = v2 / v1 
         u3 = v3 / v1 
!                                                                       
         do 210 j = km1, enorn 
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j) 
            a(k,j) = a(k,j) + t * v1 
            a(k1,j) = a(k1,j) + t * v2 
            a(k2,j) = a(k2,j) + t * v3 
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j) 
            b(k,j) = b(k,j) + t * v1 
            b(k1,j) = b(k1,j) + t * v2 
            b(k2,j) = b(k2,j) + t * v3 
  210    continue 
!                                                                       
         if (k  ==  l) go to 220 
         a(k1,km1) = 0.0_wp 
         a(k2,km1) = 0.0_wp 
!     .......... zero b(k+2,k+1) and b(k+2,k) ..........                
  220    s = abs(b(k2,k2)) + abs(b(k2,k1)) + abs(b(k2,k)) 
         if (s  ==  0.0_wp) go to 240 
         u1 = b(k2,k2) / s 
         u2 = b(k2,k1) / s 
         u3 = b(k2,k) / s 
         r = sign(sqrt(u1*u1+u2*u2+u3*u3),u1) 
         v1 = -(u1 + r) / r 
         v2 = -u2 / r 
         v3 = -u3 / r 
         u2 = v2 / v1 
         u3 = v3 / v1 
!                                                                       
         do 230 i = lor1, ll 
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k) 
            a(i,k2) = a(i,k2) + t * v1 
            a(i,k1) = a(i,k1) + t * v2 
            a(i,k) = a(i,k) + t * v3 
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k) 
            b(i,k2) = b(i,k2) + t * v1 
            b(i,k1) = b(i,k1) + t * v2 
            b(i,k) = b(i,k) + t * v3 
  230    continue 
!                                                                       
         b(k2,k) = 0.0_wp 
         b(k2,k1) = 0.0_wp 
         if (.not. matz) go to 240 
!                                                                       
         do 235 i = 1, n 
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k) 
            z(i,k2) = z(i,k2) + t * v1 
            z(i,k1) = z(i,k1) + t * v2 
            z(i,k) = z(i,k) + t * v3 
  235    continue 
!     .......... zero b(k+1,k) ..........                               
  240    s = abs(b(k1,k1)) + abs(b(k1,k)) 
         if (s  ==  0.0_wp) go to 260 
         u1 = b(k1,k1) / s 
         u2 = b(k1,k) / s 
         r = sign(sqrt(u1*u1+u2*u2),u1) 
         v1 = -(u1 + r) / r 
         v2 = -u2 / r 
         u2 = v2 / v1 
!                                                                       
         do 250 i = lor1, ll 
            t = a(i,k1) + u2 * a(i,k) 
            a(i,k1) = a(i,k1) + t * v1 
            a(i,k) = a(i,k) + t * v2 
            t = b(i,k1) + u2 * b(i,k) 
            b(i,k1) = b(i,k1) + t * v1 
            b(i,k) = b(i,k) + t * v2 
  250    continue 
!                                                                       
         b(k1,k) = 0.0_wp 
         if (.not. matz) go to 260 
!                                                                       
         do 255 i = 1, n 
            t = z(i,k1) + u2 * z(i,k) 
            z(i,k1) = z(i,k1) + t * v1 
            z(i,k) = z(i,k) + t * v2 
  255    continue 
!                                                                       
  260 continue 
!     .......... end qz step ..........                                 
      go to 70 
!     .......... set error -- all eigenvalues have not                  
!                converged after 30*n iterations ..........             
 1000 ierr = en 
!     .......... save epsb for use by qzval and qzvec ..........        
 1001 if (n  >  1) b(n,1) = epsb 
      if(epsb >= 1.0e-10_wp)                                                &
        write(*,*)'Watch the eigenvalues:  Error = ',epsb
      return 
      end subroutine qzit

!=============================================================================80

      subroutine qzval(nm,n,a,b,alfr,alfi,beta,matz,z) 
!                                                                       
      integer i,j,n,en,na,nm,nn,isw 
      real(wp) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n) 
      real(wp) c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1, &
     &       u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr,      &
     &       ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb         
      logical matz 
!                                                                       
!     this subroutine is the third step of the qz algorithm             
!     for solving generalized matrix eigenvalue problems,               
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.      
!                                                                       
!     this subroutine accepts a pair of real matrices, one of them      
!     in quasi-triangular form and the other in upper triangular form.  
!     it reduces the quasi-triangular matrix further, so that any       
!     remaining 2-by-2 blocks correspond to pairs of complex            
!     eigenvalues, and returns quantities whose ratios give the         
!     generalized eigenvalues.  it is usually preceded by  qzhes        
!     and  qzit  and may be followed by  qzvec.                         
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrices.                                
!                                                                       
!        a contains a real upper quasi-triangular matrix.               
!                                                                       
!        b contains a real upper triangular matrix.  in addition,       
!          location b(n,1) contains the tolerance quantity (epsb)       
!          computed and saved in  qzit.                                 
!                                                                       
!        matz should be set to .true. if the right hand transformations 
!          are to be accumulated for later use in computing             
!          eigenvectors, and to .false. otherwise.                      
!                                                                       
!        z contains, if matz has been set to .true., the                
!          transformation matrix produced in the reductions by qzhes    
!          and qzit, if performed, or else the identity matrix.         
!          if matz has been set to .false., z is not referenced.        
!                                                                       
!     on output                                                         
!                                                                       
!        a has been reduced further to a quasi-triangular matrix        
!          in which all nonzero subdiagonal elements correspond to      
!          pairs of complex eigenvalues.                                
!                                                                       
!        b is still in upper triangular form, although its elements     
!          have been altered.  b(n,1) is unaltered.                     
!                                                                       
!        alfr and alfi contain the real and imaginary parts of the      
!          diagonal elements of the triangular matrix that would be     
!          obtained if a were reduced completely to triangular form     
!          by unitary transformations.  non-zero values of alfi occur   
!          in pairs, the first member positive and the second negative. 
!                                                                       
!        beta contains the diagonal elements of the corresponding b,    
!          normalized to be real and non-negative.  the generalized     
!          eigenvalues are then the ratios ((alfr+i*alfi)/beta).        
!                                                                       
!        z contains the product of the right hand transformations       
!          (for all three steps) if matz has been set to .true.         
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      epsb = b(n,1) 
      isw = 1 
!     .......... find eigenvalues of quasi-triangular matrices.         
!                for en=n step -1 until 1 do -- ..........              
      do 510 nn = 1, n 
         en = n + 1 - nn 
         na = en - 1 
         if (isw  ==  2) go to 505 
         if (en  ==  1) go to 410 
         if (a(en,na)  /=  0.0_wp) go to 420 
!     .......... 1-by-1 block, one real root ..........                 
  410    alfr(en) = a(en,en) 
         if (b(en,en)  <  0.0_wp) alfr(en) = -alfr(en) 
         beta(en) = abs(b(en,en)) 
         alfi(en) = 0.0_wp 
         go to 510 
!     .......... 2-by-2 block ..........                                
  420    if (abs(b(na,na))  <=  epsb) go to 455 
         if (abs(b(en,en))  >  epsb) go to 430 
         a1 = a(en,en) 
         a2 = a(en,na) 
         bn = 0.0_wp 
         go to 435 
  430    an = abs(a(na,na)) + abs(a(na,en)) + abs(a(en,na))          &
     &      + abs(a(en,en))                                            
         bn = abs(b(na,na)) + abs(b(na,en)) + abs(b(en,en)) 
         a11 = a(na,na) / an 
         a12 = a(na,en) / an 
         a21 = a(en,na) / an 
         a22 = a(en,en) / an 
         b11 = b(na,na) / bn 
         b12 = b(na,en) / bn 
         b22 = b(en,en) / bn 
         e = a11 / b11 
         ei = a22 / b22 
         s = a21 / (b11 * b22) 
         t = (a22 - e * b22) / b22 
         if (abs(e)  <=  abs(ei)) go to 431 
         e = ei 
         t = (a11 - e * b11) / b11 
  431    c = 0.5_wp * (t - s * b12) 
         d = c * c + s * (a12 - e * b12) 
         if (d  <  0.0_wp) go to 480 
!     .......... two real roots.                                        
!                zero both a(en,na) and b(en,na) ..........             
         e = e + (c + sign(sqrt(d),c)) 
         a11 = a11 - e * b11 
         a12 = a12 - e * b12 
         a22 = a22 - e * b22 
         if (abs(a11) + abs(a12)  <                                  &
     &       abs(a21) + abs(a22)) go to 432                           
         a1 = a12 
         a2 = a11 
         go to 435 
  432    a1 = a22 
         a2 = a21 
!     .......... choose and apply real z ..........                     
  435    s = abs(a1) + abs(a2) 
         u1 = a1 / s 
         u2 = a2 / s 
         r = sign(sqrt(u1*u1+u2*u2),u1) 
         v1 = -(u1 + r) / r 
         v2 = -u2 / r 
         u2 = v2 / v1 
!                                                                       
         do 440 i = 1, en 
            t = a(i,en) + u2 * a(i,na) 
            a(i,en) = a(i,en) + t * v1 
            a(i,na) = a(i,na) + t * v2 
            t = b(i,en) + u2 * b(i,na) 
            b(i,en) = b(i,en) + t * v1 
            b(i,na) = b(i,na) + t * v2 
  440    continue 
!                                                                       
         if (.not. matz) go to 450 
!                                                                       
         do 445 i = 1, n 
            t = z(i,en) + u2 * z(i,na) 
            z(i,en) = z(i,en) + t * v1 
            z(i,na) = z(i,na) + t * v2 
  445    continue 
!                                                                       
  450    if (bn  ==  0.0_wp) go to 475 
         if (an  <  abs(e) * bn) go to 455 
         a1 = b(na,na) 
         a2 = b(en,na) 
         go to 460 
  455    a1 = a(na,na) 
         a2 = a(en,na) 
!     .......... choose and apply real q ..........                     
  460    s = abs(a1) + abs(a2) 
         if (s  ==  0.0_wp) go to 475 
         u1 = a1 / s 
         u2 = a2 / s 
         r = sign(sqrt(u1*u1+u2*u2),u1) 
         v1 = -(u1 + r) / r 
         v2 = -u2 / r 
         u2 = v2 / v1 
!                                                                       
         do 470 j = na, n 
            t = a(na,j) + u2 * a(en,j) 
            a(na,j) = a(na,j) + t * v1 
            a(en,j) = a(en,j) + t * v2 
            t = b(na,j) + u2 * b(en,j) 
            b(na,j) = b(na,j) + t * v1 
            b(en,j) = b(en,j) + t * v2 
  470    continue 
!                                                                       
  475    a(en,na) = 0.0_wp 
         b(en,na) = 0.0_wp 
         alfr(na) = a(na,na) 
         alfr(en) = a(en,en) 
         if (b(na,na)  <  0.0_wp) alfr(na) = -alfr(na) 
         if (b(en,en)  <  0.0_wp) alfr(en) = -alfr(en) 
         beta(na) = abs(b(na,na)) 
         beta(en) = abs(b(en,en)) 
         alfi(en) = 0.0_wp 
         alfi(na) = 0.0_wp 
         go to 505 
!     .......... two complex roots ..........                           
  480    e = e + c 
         ei = sqrt(-d) 
         a11r = a11 - e * b11 
         a11i = ei * b11 
         a12r = a12 - e * b12 
         a12i = ei * b12 
         a22r = a22 - e * b22 
         a22i = ei * b22 
         if (abs(a11r) + abs(a11i) + abs(a12r) + abs(a12i)  <      &
     &       abs(a21) + abs(a22r) + abs(a22i)) go to 482             
         a1 = a12r 
         a1i = a12i 
         a2 = -a11r 
         a2i = -a11i 
         go to 485 
  482    a1 = a22r 
         a1i = a22i 
         a2 = -a21 
         a2i = 0.0_wp 
!     .......... choose complex z ..........                            
  485    cz = sqrt(a1*a1+a1i*a1i) 
         if (cz  ==  0.0_wp) go to 487 
         szr = (a1 * a2 + a1i * a2i) / cz 
         szi = (a1 * a2i - a1i * a2) / cz 
         r = sqrt(cz*cz+szr*szr+szi*szi) 
         cz = cz / r 
         szr = szr / r 
         szi = szi / r 
         go to 490 
  487    szr = 1.0_wp 
         szi = 0.0_wp 
  490    if (an  <  (abs(e) + ei) * bn) go to 492 
         a1 = cz * b11 + szr * b12 
         a1i = szi * b12 
         a2 = szr * b22 
         a2i = szi * b22 
         go to 495 
  492    a1 = cz * a11 + szr * a12 
         a1i = szi * a12 
         a2 = cz * a21 + szr * a22 
         a2i = szi * a22 
!     .......... choose complex q ..........                            
  495    cq = sqrt(a1*a1+a1i*a1i) 
         if (cq  ==  0.0_wp) go to 497 
         sqr = (a1 * a2 + a1i * a2i) / cq 
         sqi = (a1 * a2i - a1i * a2) / cq 
         r = sqrt(cq*cq+sqr*sqr+sqi*sqi) 
         cq = cq / r 
         sqr = sqr / r 
         sqi = sqi / r 
         go to 500 
  497    sqr = 1.0_wp 
         sqi = 0.0_wp 
!     .......... compute diagonal elements that would result            
!                if transformations were applied ..........             
  500    ssr = sqr * szr + sqi * szi 
         ssi = sqr * szi - sqi * szr 
         i = 1 
         tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21           &
     &      + ssr * a22                                                 
         ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22 
         dr = cq * cz * b11 + cq * szr * b12 + ssr * b22 
         di = cq * szi * b12 + ssi * b22 
         go to 503 
  502    i = 2 
         tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21               &
     &      + cq * cz * a22                                             
         ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21 
         dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22 
         di = -ssi * b11 - sqi * cz * b12 
  503    t = ti * dr - tr * di 
         j = na 
         if (t  <  0.0_wp) j = en 
         r = sqrt(dr*dr+di*di) 
         beta(j) = bn * r 
         alfr(j) = an * (tr * dr + ti * di) / r 
         alfi(j) = an * t / r 
         if (i  ==  1) go to 502 
  505    isw = 3 - isw 
  510 continue 
      b(n,1) = epsb 
!                                                                       
      return 
      end subroutine qzval

!=============================================================================80

      subroutine qzvec(nm,n,a,b,alfr,alfi,beta,z) 
!                                                                       
      integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2 
      real(wp) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n) 
      real(wp) d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1,&
     &       zz,z1,alfm,almi,almr,betm,epsb                             
!                                                                       
!     this subroutine is the optional fourth step of the qz algorithm   
!     for solving generalized matrix eigenvalue problems,               
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.      
!                                                                       
!     this subroutine accepts a pair of real matrices, one of them in   
!     quasi-triangular form (in which each 2-by-2 block corresponds to  
!     a pair of complex eigenvalues) and the other in upper triangular  
!     form.  it computes the eigenvectors of the triangular problem and 
!     transforms the results back to the original coordinate system.    
!     it is usually preceded by  qzhes,  qzit, and  qzval.              
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrices.                                
!                                                                       
!        a contains a real upper quasi-triangular matrix.               
!                                                                       
!        b contains a real upper triangular matrix.  in addition,       
!          location b(n,1) contains the tolerance quantity (epsb)       
!          computed and saved in  qzit.                                 
!                                                                       
!        alfr, alfi, and beta  are vectors with components whose        
!          ratios ((alfr+i*alfi)/beta) are the generalized              
!          eigenvalues.  they are usually obtained from  qzval.         
!                                                                       
!        z contains the transformation matrix produced in the           
!          reductions by  qzhes,  qzit, and  qzval, if performed.       
!          if the eigenvectors of the triangular problem are            
!          desired, z must contain the identity matrix.                 
!                                                                       
!     on output                                                         
!                                                                       
!        a is unaltered.  its subdiagonal elements provide information  
!           about the storage of the complex eigenvectors.              
!                                                                       
!        b has been destroyed.                                          
!                                                                       
!        alfr, alfi, and beta are unaltered.                            
!                                                                       
!        z contains the real and imaginary parts of the eigenvectors.   
!          if alfi(i)  ==  0.0, the i-th eigenvalue is real and         
!            the i-th column of z contains its eigenvector.             
!          if alfi(i)  /=  0.0, the i-th eigenvalue is complex.         
!            if alfi(i)  >  0.0, the eigenvalue is the first of        
!              a complex pair and the i-th and (i+1)-th columns         
!              of z contain its eigenvector.                            
!            if alfi(i)  <  0.0, the eigenvalue is the second of       
!              a complex pair and the (i-1)-th and i-th columns         
!              of z contain the conjugate of its eigenvector.           
!          each eigenvector is normalized so that the modulus           
!          of its largest component is 1.0 .                            
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      epsb = b(n,1) 
      isw = 1 
!     .......... for en=n step -1 until 1 do -- ..........              
      do 800 nn = 1, n 
         en = n + 1 - nn 
         na = en - 1 
         if (isw  ==  2) go to 795 
         if (alfi(en)  /=  0.0_wp) go to 710 
!     .......... real vector ..........                                 
         m = en 
         b(en,en) = 1.0_wp 
         if (na  ==  0) go to 800 
         alfm = alfr(m) 
         betm = beta(m) 
!     .......... for i=en-1 step -1 until 1 do -- ..........            
         do 700 ii = 1, na 
            i = en - ii 
            w = betm * a(i,i) - alfm * b(i,i) 
            r = 0.0_wp 
!                                                                       
            do 610 j = m, en 
  610       r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en) 
!                                                                       
            if (i  ==  1 .or. isw  ==  2) go to 630 
            if (betm * a(i,i-1)  ==  0.0_wp) go to 630 
            zz = w 
            s = r 
            go to 690 
  630       m = i 
            if (isw  ==  2) go to 640 
!     .......... real 1-by-1 block ..........                           
            t = w 
            if (w  ==  0.0_wp) t = epsb 
            b(i,en) = -r / t 
            go to 700 
!     .......... real 2-by-2 block ..........                           
  640       x = betm * a(i,i+1) - alfm * b(i,i+1) 
            y = betm * a(i+1,i) 
            q = w * zz - x * y 
            t = (x * s - zz * r) / q 
            b(i,en) = t 
            if (abs(x)  <=  abs(zz)) go to 650 
            b(i+1,en) = (-r - w * t) / x 
            go to 690 
  650       b(i+1,en) = (-s - y * t) / zz 
  690       isw = 3 - isw 
  700    continue 
!     .......... end real vector ..........                             
         go to 800 
!     .......... complex vector ..........                              
  710    m = na 
         almr = alfr(m) 
         almi = alfi(m) 
         betm = beta(m) 
!     .......... last vector component chosen imaginary so that         
!                eigenvector matrix is triangular ..........            
         y = betm * a(en,na) 
         b(na,na) = -almi * b(en,en) / y 
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y 
         b(en,na) = 0.0_wp 
         b(en,en) = 1.0_wp 
         enm2 = na - 1 
         if (enm2  ==  0) go to 795 
!     .......... for i=en-2 step -1 until 1 do -- ..........            
         do 790 ii = 1, enm2 
            i = na - ii 
            w = betm * a(i,i) - almr * b(i,i) 
            w1 = -almi * b(i,i) 
            ra = 0.0_wp 
            sa = 0.0_wp 
!                                                                       
            do 760 j = m, en 
               x = betm * a(i,j) - almr * b(i,j) 
               x1 = -almi * b(i,j) 
               ra = ra + x * b(j,na) - x1 * b(j,en) 
               sa = sa + x * b(j,en) + x1 * b(j,na) 
  760       continue 
!                                                                       
            if (i  ==  1 .or. isw  ==  2) go to 770 
            if (betm * a(i,i-1)  ==  0.0_wp) go to 770 
            zz = w 
            z1 = w1 
            r = ra 
            s = sa 
            isw = 2 
            go to 790 
  770       m = i 
            if (isw  ==  2) go to 780 
!     .......... complex 1-by-1 block ..........                        
            tr = -ra 
            ti = -sa 
  773       dr = w 
            di = w1 
!     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........  
  775       if (abs(di)  >  abs(dr)) go to 777 
            rr = di / dr 
            d = dr + di * rr 
            t1 = (tr + ti * rr) / d 
            t2 = (ti - tr * rr) / d 
            go to (787,782), isw 
  777       rr = dr / di 
            d = dr * rr + di 
            t1 = (tr * rr + ti) / d 
            t2 = (ti * rr - tr) / d 
            go to (787,782), isw 
!     .......... complex 2-by-2 block ..........                        
  780       x = betm * a(i,i+1) - almr * b(i,i+1) 
            x1 = -almi * b(i,i+1) 
            y = betm * a(i+1,i) 
            tr = y * ra - w * r + w1 * s 
            ti = y * sa - w * s - w1 * r 
            dr = w * zz - w1 * z1 - x * y 
            di = w * z1 + w1 * zz - x1 * y 
            if (dr  ==  0.0_wp .and. di  ==  0.0_wp) dr = epsb 
            go to 775 
  782       b(i+1,na) = t1 
            b(i+1,en) = t2 
            isw = 1 
            if (abs(y)  >  abs(w) + abs(w1)) go to 785 
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en) 
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na) 
            go to 773 
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y 
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y 
  787       b(i,na) = t1 
            b(i,en) = t2 
  790    continue 
!     .......... end complex vector ..........                          
  795    isw = 3 - isw 
  800 continue 
!     .......... end back substitution.                                 
!                transform to original coordinate system.               
!                for j=n step -1 until 1 do -- ..........               
      do 880 jj = 1, n 
         j = n + 1 - jj 
!                                                                       
         do 880 i = 1, n 
            zz = 0.0_wp 
!                                                                       
            do 860 k = 1, j 
  860       zz = zz + z(i,k) * b(k,j) 
!                                                                       
            z(i,j) = zz 
  880 continue 
!     .......... normalize so that modulus of largest                   
!                component of each vector is 1.                         
!                (isw is 1 initially from before) ..........            
      do 950 j = 1, n 
         d = 0.0_wp 
         if (isw  ==  2) go to 920 
         if (alfi(j)  /=  0.0_wp) go to 945 
!                                                                       
         do 890 i = 1, n 
            if (abs(z(i,j))  >  d) d = abs(z(i,j)) 
  890    continue 
!                                                                       
         do 900 i = 1, n 
  900    z(i,j) = z(i,j) / d 
!                                                                       
         go to 950 
!                                                                       
  920    do 930 i = 1, n 
            r = abs(z(i,j-1)) + abs(z(i,j)) 
            if (r  /=  0.0_wp) r = r * sqrt((z(i,j-1)/r)**2             &
     &                                     +(z(i,j)/r)**2)              
            if (r  >  d) d = r 
  930    continue 
!                                                                       
         do 940 i = 1, n 
            z(i,j-1) = z(i,j-1) / d 
            z(i,j) = z(i,j) / d 
  940    continue 
!                                                                       
  945    isw = 3 - isw 
  950 continue 
!                                                                       
      return 
      end subroutine qzvec

!=============================================================================80

      subroutine rebak(nm,n,b,dl,m,z) 
!                                                                       
      integer i,j,k,m,n,i1,ii,nm 
      real(wp) b(nm,n),dl(n),z(nm,m) 
      real(wp) x 
!                                                                       
!     this subroutine is a translation of the algol procedure rebaka,   
!     num. math. 11, 99-110(1968) by martin and wilkinson.              
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).   
!                                                                       
!     this subroutine forms the eigenvectors of a generalized           
!     symmetric eigensystem by back transforming those of the           
!     derived symmetric matrix determined by  reduc.                    
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix system.                           
!                                                                       
!        b contains information about the similarity transformation     
!          (cholesky decomposition) used in the reduction by  reduc     
!          in its strict lower triangle.                                
!                                                                       
!        dl contains further information about the transformation.      
!                                                                       
!        m is the number of eigenvectors to be back transformed.        
!                                                                       
!        z contains the eigenvectors to be back transformed             
!          in its first m columns.                                      
!                                                                       
!     on output                                                         
!                                                                       
!        z contains the transformed eigenvectors                        
!          in its first m columns.                                      
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (m  ==  0) go to 200 
!                                                                       
      do 100 j = 1, m 
!     .......... for i=n step -1 until 1 do -- ..........               
         do 100 ii = 1, n 
            i = n + 1 - ii 
            i1 = i + 1 
            x = z(i,j) 
            if (i  ==  n) go to 80 
!                                                                       
            do 60 k = i1, n 
   60       x = x - b(k,i) * z(k,j) 
!                                                                       
   80       z(i,j) = x / dl(i) 
  100 continue 
!                                                                       
  200 return 
      end subroutine rebak

!=============================================================================80

      subroutine reduc(nm,n,a,b,dl,ierr) 
!                                                                       
      integer i,j,k,n,i1,j1,nm,nn,ierr 
      real(wp) a(nm,n),b(nm,n),dl(n) 
      real(wp) x,y 
!                                                                       
!     this subroutine is a translation of the algol procedure reduc1,   
!     num. math. 11, 99-110(1968) by martin and wilkinson.              
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).   
!                                                                       
!     this subroutine reduces the generalized symmetric eigenproblem    
!     ax=(lambda)bx, where b is positive definite, to the standard      
!     symmetric eigenproblem using the cholesky factorization of b.     
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrices a and b.  if the cholesky       
!          factor l of b is already available, n should be prefixed     
!          with a minus sign.                                           
!                                                                       
!        a and b contain the real symmetric input matrices.  only the   
!          full upper triangles of the matrices need be supplied.  if   
!          n is negative, the strict lower triangle of b contains,      
!          instead, the strict lower triangle of its cholesky factor l. 
!                                                                       
!        dl contains, if n is negative, the diagonal elements of l.     
!                                                                       
!     on output                                                         
!                                                                       
!        a contains in its full lower triangle the full lower triangle  
!          of the symmetric matrix derived from the reduction to the    
!          standard form.  the strict upper triangle of a is unaltered. 
!                                                                       
!        b contains in its strict lower triangle the strict lower       
!          triangle of its cholesky factor l.  the full upper           
!          triangle of b is unaltered.                                  
!                                                                       
!        dl contains the diagonal elements of l.                        
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          7*n+1      if b is not positive definite.                    
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      nn = iabs(n) 
      if (n  <  0) go to 100 
!     .......... form l in the arrays b and dl ..........               
      do 80 i = 1, n 
         i1 = i - 1 
!                                                                       
         do 80 j = i, n 
            x = b(i,j) 
            if (i  ==  1) go to 40 
!                                                                       
            do 20 k = 1, i1 
   20       x = x - b(i,k) * b(j,k) 
!                                                                       
   40       if (j  /=  i) go to 60 
            if (x  <=  0.0_wp) go to 1000 
            y = sqrt(x) 
            dl(i) = y 
            go to 80 
   60       b(j,i) = x / y 
   80 continue 
!     .......... form the transpose of the upper triangle of inv(l)*a   
!                in the lower triangle of the array a ..........        
  100 do 200 i = 1, nn 
         i1 = i - 1 
         y = dl(i) 
!                                                                       
         do 200 j = i, nn 
            x = a(i,j) 
            if (i  ==  1) go to 180 
!                                                                       
            do 160 k = 1, i1 
  160       x = x - b(i,k) * a(j,k) 
!                                                                       
  180       a(j,i) = x / y 
  200 continue 
!     .......... pre-multiply by inv(l) and overwrite ..........        
      do 300 j = 1, nn 
         j1 = j - 1 
!                                                                       
         do 300 i = j, nn 
            x = a(i,j) 
            if (i  ==  j) go to 240 
            i1 = i - 1 
!                                                                       
            do 220 k = j, i1 
  220       x = x - a(k,j) * b(i,k) 
!                                                                       
  240       if (j  ==  1) go to 280 
!                                                                       
            do 260 k = 1, j1 
  260       x = x - a(j,k) * b(i,k) 
!                                                                       
  280       a(i,j) = x / dl(i) 
  300 continue 
!                                                                       
      go to 1001 
!     .......... set error -- b is not positive definite ..........     
 1000 ierr = 7 * n + 1 
 1001 return 
      end subroutine reduc

!=============================================================================80

      subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr) 
!                                                                       
      integer n,nm,is1,is2,ierr,matz 
      real(wp) a(nm,n),wr(n),wi(n),z(nm,n),fv1(n) 
      integer iv1(n) 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find the eigenvalues and eigenvectors (if desired)             
!     of a real general matrix.                                         
!                                                                       
!     on input                                                          
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement.                                           
!                                                                       
!        n  is the order of the matrix  a.                              
!                                                                       
!        a  contains the real general matrix.                           
!                                                                       
!        matz  is an integer variable set equal to zero if              
!        only eigenvalues are desired.  otherwise it is set to          
!        any non-zero integer for both eigenvalues and eigenvectors.    
!                                                                       
!     on output                                                         
!                                                                       
!        wr  and  wi  contain the real and imaginary parts,             
!        respectively, of the eigenvalues.  complex conjugate           
!        pairs of eigenvalues appear consecutively with the             
!        eigenvalue having the positive imaginary part first.           
!                                                                       
!        z  contains the real and imaginary parts of the eigenvectors   
!        if matz is not zero.  if the j-th eigenvalue is real, the      
!        j-th column of  z  contains its eigenvector.  if the j-th      
!        eigenvalue is complex with positive imaginary part, the        
!        j-th and (j+1)-th columns of  z  contain the real and          
!        imaginary parts of its eigenvector.  the conjugate of this     
!        vector is the eigenvector for the conjugate eigenvalue.        
!                                                                       
!        ierr  is an integer output variable set equal to an error      
!           completion code described in the documentation for hqr      
!           and hqr2.  the normal completion code is zero.              
!                                                                       
!        iv1  and  fv1  are temporary storage arrays.                   
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (n  <=  nm) go to 10 
      ierr = 10 * n 
      go to 50 
!                                                                       
   10 call  balanc(nm,n,a,is1,is2,fv1) 
      call  elmhes(nm,n,is1,is2,a,iv1) 
      if (matz  /=  0) go to 20 
!     .......... find eigenvalues only ..........                       
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr) 
      go to 50 
!     .......... find both eigenvalues and eigenvectors ..........      
   20 call  eltran(nm,n,is1,is2,a,iv1,z) 
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr) 
      if (ierr  /=  0) go to 50 
      call  balbak(nm,n,is1,is2,fv1,n,z) 
   50 return 
      end subroutine rg

!=============================================================================80

      subroutine rgg(nm,n,a,b,alfr,alfi,beta,matz,z,ierr) 
!                                                                       
      integer, parameter               :: dp = 8    ! double precision

      integer n,nm,ierr,matz 
      real(wp) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n) 
      logical tf 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find the eigenvalues and eigenvectors (if desired)             
!     for the real general generalized eigenproblem  ax = (lambda)bx.   
!                                                                       
!     on input                                                          
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement.                                           
!                                                                       
!        n  is the order of the matrices  a  and  b.                    
!                                                                       
!        a  contains a real general matrix.                             
!                                                                       
!        b  contains a real general matrix.                             
!                                                                       
!        matz  is an integer variable set equal to zero if              
!        only eigenvalues are desired.  otherwise it is set to          
!        any non-zero integer for both eigenvalues and eigenvectors.    
!                                                                       
!     on output                                                         
!                                                                       
!        alfr  and  alfi  contain the real and imaginary parts,         
!        respectively, of the numerators of the eigenvalues.            
!                                                                       
!        beta  contains the denominators of the eigenvalues,            
!        which are thus given by the ratios  (alfr+i*alfi)/beta.        
!        complex conjugate pairs of eigenvalues appear consecutively    
!        with the eigenvalue having the positive imaginary part first.  
!                                                                       
!        z  contains the real and imaginary parts of the eigenvectors   
!        if matz is not zero.  if the j-th eigenvalue is real, the      
!        j-th column of  z  contains its eigenvector.  if the j-th      
!        eigenvalue is complex with positive imaginary part, the        
!        j-th and (j+1)-th columns of  z  contain the real and          
!        imaginary parts of its eigenvector.  the conjugate of this     
!        vector is the eigenvector for the conjugate eigenvalue.        
!                                                                       
!        ierr  is an integer output variable set equal to an error      
!           completion code described in the documentation for qzit.    
!           the normal completion code is zero.                         
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (n  <=  nm) go to 10 
      ierr = 10 * n 
      go to 50 
!                                                                       
   10 if (matz  /=  0) go to 20 
!     .......... find eigenvalues only ..........                       
      tf = .false. 
      call  qzhes(nm,n,a,b,tf,z) 
      call  qzit(nm,n,a,b,0.0_wp,tf,z,ierr) 
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z) 
      go to 50 
!     .......... find both eigenvalues and eigenvectors ..........      
   20 tf = .true. 
      call  qzhes(nm,n,a,b,tf,z) 
      call  qzit(nm,n,a,b,0.0_wp,tf,z,ierr) 
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z) 
      if (ierr  /=  0) go to 50 
      call  qzvec(nm,n,a,b,alfr,alfi,beta,z) 
   50 return 
      end subroutine rgg

!=============================================================================80

      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr) 
!                                                                       
      integer n,nm,ierr,matz 
      real(wp) a(nm,n),w(n),z(nm,n),fv1(n),fv2(n) 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find the eigenvalues and eigenvectors (if desired)             
!     of a real symmetric matrix.                                       
!                                                                       
!     on input                                                          
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement.                                           
!                                                                       
!        n  is the order of the matrix  a.                              
!                                                                       
!        a  contains the real symmetric matrix.                         
!                                                                       
!        matz  is an integer variable set equal to zero if              
!        only eigenvalues are desired.  otherwise it is set to          
!        any non-zero integer for both eigenvalues and eigenvectors.    
!                                                                       
!     on output                                                         
!                                                                       
!        w  contains the eigenvalues in ascending order.                
!                                                                       
!        z  contains the eigenvectors if matz is not zero.              
!                                                                       
!        ierr  is an integer output variable set equal to an error      
!           completion code described in the documentation for tqlrat   
!           and tql2.  the normal completion code is zero.              
!                                                                       
!        fv1  and  fv2  are temporary storage arrays.                   
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (n  <=  nm) go to 10 
      ierr = 10 * n 
      go to 50 
!                                                                       
   10 if (matz  /=  0) go to 20 
!     .......... find eigenvalues only ..........                       
      call  tred1(nm,n,a,w,fv1,fv2) 
!  tqlrat encounters catastrophic underflow on the Vax                  
!     call  tqlrat(n,w,fv2,ierr)                                        
      call  tql1(n,w,fv1,ierr) 
      go to 50 
!     .......... find both eigenvalues and eigenvectors ..........      
   20 call  tred2(nm,n,a,w,fv1,z) 
      call  tql2(nm,n,w,fv1,z,ierr) 
   50 return 
      end subroutine rs

!=============================================================================80

      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr) 
!                                                                       
      integer n,nm,ierr,matz 
      real(wp) a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n) 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find the eigenvalues and eigenvectors (if desired)             
!     for the real symmetric generalized eigenproblem  ax = (lambda)bx. 
!                                                                       
!     on input                                                          
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement.                                           
!                                                                       
!        n  is the order of the matrices  a  and  b.                    
!                                                                       
!        a  contains a real symmetric matrix.                           
!                                                                       
!        b  contains a positive definite real symmetric matrix.         
!                                                                       
!        matz  is an integer variable set equal to zero if              
!        only eigenvalues are desired.  otherwise it is set to          
!        any non-zero integer for both eigenvalues and eigenvectors.    
!                                                                       
!     on output                                                         
!                                                                       
!        w  contains the eigenvalues in ascending order.                
!                                                                       
!        z  contains the eigenvectors if matz is not zero.              
!                                                                       
!        ierr  is an integer output variable set equal to an error      
!           completion code described in the documentation for tqlrat   
!           and tql2.  the normal completion code is zero.              
!                                                                       
!        fv1  and  fv2  are temporary storage arrays.                   
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (n  <=  nm) go to 10 
      ierr = 10 * n 
      go to 50 
!                                                                       
   10 call  reduc(nm,n,a,b,fv2,ierr) 
      if (ierr  /=  0) go to 50 
      if (matz  /=  0) go to 20 
!     .......... find eigenvalues only ..........                       
      call  tred1(nm,n,a,w,fv1,fv2) 
      call  tqlrat(n,w,fv2,ierr) 
      go to 50 
!     .......... find both eigenvalues and eigenvectors ..........      
   20 call  tred2(nm,n,a,w,fv1,z) 
      call  tql2(nm,n,w,fv1,z,ierr) 
      if (ierr  /=  0) go to 50 
      call  rebak(nm,n,b,fv2,n,z) 
   50 return 
      end subroutine rsg

!=============================================================================80

      subroutine rsm(nm,n,a,w,m,z,fwork,iwork,ierr) 
!                                                                       
      integer n,nm,m,iwork(n),ierr 
      integer k1,k2,k3,k4,k5,k6,k7,k8
      real(wp) a(nm,n),w(n),z(nm,m),fwork(1) 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find all of the eigenvalues and some of the eigenvectors       
!     of a real symmetric matrix.                                       
!                                                                       
!     on input                                                          
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement.                                           
!                                                                       
!        n  is the order of the matrix  a.                              
!                                                                       
!        a  contains the real symmetric matrix.                         
!                                                                       
!        m  the eigenvectors corresponding to the first m eigenvalues   
!           are to be computed.                                         
!           if m = 0 then no eigenvectors are computed.                 
!           if m = n then all of the eigenvectors are computed.         
!                                                                       
!     on output                                                         
!                                                                       
!        w  contains all n eigenvalues in ascending order.              
!                                                                       
!        z  contains the orthonormal eigenvectors associated with       
!           the first m eigenvalues.                                    
!                                                                       
!        ierr  is an integer output variable set equal to an error      
!           completion code described in the documentation for tqlrat,  
!           imtqlv and tinvit.  the normal completion code is zero.     
!                                                                       
!        fwork  is a temporary storage array of dimension 8*n.          
!                                                                       
!        iwork  is an integer temporary storage array of dimension n.   
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 10 * n 
      if (n  >  nm .or. m  >  nm) go to 50 
      k1 = 1 
      k2 = k1 + n 
      k3 = k2 + n 
      k4 = k3 + n 
      k5 = k4 + n 
      k6 = k5 + n 
      k7 = k6 + n 
      k8 = k7 + n 
      if (m  >  0) go to 10 
!     .......... find eigenvalues only ..........                       
      call  tred1(nm,n,a,w,fwork(k1),fwork(k2)) 
      call  tqlrat(n,w,fwork(k2),ierr) 
      go to 50 
!     .......... find all eigenvalues and m eigenvectors ..........     
   10 call  tred1(nm,n,a,fwork(k1),fwork(k2),fwork(k3)) 
      call  imtqlv(n,fwork(k1),fwork(k2),fwork(k3),w,iwork,             &
     &             ierr,fwork(k4))                                      
      call  tinvit(nm,n,fwork(k1),fwork(k2),fwork(k3),m,w,iwork,z,ierr, &
     &             fwork(k4),fwork(k5),fwork(k6),fwork(k7),fwork(k8))   
      call  trbak1(nm,n,a,fwork(k2),m,z) 
   50 return 
      end subroutine rsm

!=============================================================================80

      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1) 
!                                                                       
      integer, parameter               :: dp = 8    ! double precision

      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr 
      real(wp) a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n) 
      real(wp) c,f,g,h,s,x,y,z,tst1,tst2,scale
      logical matu,matv 
!                                                                       
!     this subroutine is a translation of the algol procedure svd,      
!     num. math. 14, 403-420(1970) by golub and reinsch.                
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).   
!                                                                       
!     this subroutine determines the singular value decomposition       
!          t                                                            
!     a=usv  of a real m by n rectangular matrix.  householder          
!     bidiagonalization and a variant of the qr algorithm are used.     
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.  note that nm must be at least          
!          as large as the maximum of m and n.                          
!                                                                       
!        m is the number of rows of a (and u).                          
!                                                                       
!        n is the number of columns of a (and u) and the order of v.    
!                                                                       
!        a contains the rectangular input matrix to be decomposed.      
!                                                                       
!        matu should be set to .true. if the u matrix in the            
!          decomposition is desired, and to .false. otherwise.          
!                                                                       
!        matv should be set to .true. if the v matrix in the            
!          decomposition is desired, and to .false. otherwise.          
!                                                                       
!     on output                                                         
!                                                                       
!        a is unaltered (unless overwritten by u or v).                 
!                                                                       
!        w contains the n (non-negative) singular values of a (the      
!          diagonal elements of s).  they are unordered.  if an         
!          error exit is made, the singular values should be correct    
!          for indices ierr+1,ierr+2,...,n.                             
!                                                                       
!        u contains the matrix u (orthogonal column vectors) of the     
!          decomposition if matu has been set to .true.  otherwise      
!          u is used as a temporary array.  u may coincide with a.      
!          if an error exit is made, the columns of u corresponding     
!          to indices of correct singular values should be correct.     
!                                                                       
!        v contains the matrix v (orthogonal) of the decomposition if   
!          matv has been set to .true.  otherwise v is not referenced.  
!          v may also coincide with a if u is not needed.  if an error  
!          exit is made, the columns of v corresponding to indices of   
!          correct singular values should be correct.                   
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          k          if the k-th singular value has not been           
!                     determined after 30 iterations.                   
!                                                                       
!        rv1 is a temporary storage array.                              
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
!                                                                       
      do 100 i = 1, m 
!                                                                       
         do 100 j = 1, n 
            u(i,j) = a(i,j) 
  100 continue 
!     .......... householder reduction to bidiagonal form ..........    
      g = 0.0_wp 
      scale = 0.0_wp 
      x = 0.0_wp 
!                                                                       
      do 300 i = 1, n 
         l = i + 1 
         rv1(i) = scale * g 
         g = 0.0_wp 
         s = 0.0_wp 
         scale = 0.0_wp 
         if (i  >  m) go to 210 
!                                                                       
         do 120 k = i, m 
  120    scale = scale + abs(u(k,i)) 
!                                                                       
         if (scale  ==  0.0_wp) go to 210 
!                                                                       
         do 130 k = i, m 
            u(k,i) = u(k,i) / scale 
            s = s + u(k,i)**2 
  130    continue 
!                                                                       
         f = u(i,i) 
         g = -sign(sqrt(s),f) 
         h = f * g - s 
         u(i,i) = f - g 
         if (i  ==  n) go to 190 
!                                                                       
         do 150 j = l, n 
            s = 0.0_wp 
!                                                                       
            do 140 k = i, m 
  140       s = s + u(k,i) * u(k,j) 
!                                                                       
            f = s / h 
!                                                                       
            do 150 k = i, m 
               u(k,j) = u(k,j) + f * u(k,i) 
  150    continue 
!                                                                       
  190    do 200 k = i, m 
  200    u(k,i) = scale * u(k,i) 
!                                                                       
  210    w(i) = scale * g 
         g = 0.0_wp 
         s = 0.0_wp 
         scale = 0.0_wp 
         if (i  >  m .or. i  ==  n) go to 290 
!                                                                       
         do 220 k = l, n 
  220    scale = scale + abs(u(i,k)) 
!                                                                       
         if (scale  ==  0.0_wp) go to 290 
!                                                                       
         do 230 k = l, n 
            u(i,k) = u(i,k) / scale 
            s = s + u(i,k)**2 
  230    continue 
!                                                                       
         f = u(i,l) 
         g = -sign(sqrt(s),f) 
         h = f * g - s 
         u(i,l) = f - g 
!                                                                       
         do 240 k = l, n 
  240    rv1(k) = u(i,k) / h 
!                                                                       
         if (i  ==  m) go to 270 
!                                                                       
         do 260 j = l, m 
            s = 0.0_wp 
!                                                                       
            do 250 k = l, n 
  250       s = s + u(j,k) * u(i,k) 
!                                                                       
            do 260 k = l, n 
               u(j,k) = u(j,k) + s * rv1(k) 
  260    continue 
!                                                                       
  270    do 280 k = l, n 
  280    u(i,k) = scale * u(i,k) 
!                                                                       
  290    x = max(x,abs(w(i))+abs(rv1(i))) 
  300 continue 
!     .......... accumulation of right-hand transformations ..........  
      if (.not. matv) go to 410 
!     .......... for i=n step -1 until 1 do -- ..........               
      do 400 ii = 1, n 
         i = n + 1 - ii 
         if (i  ==  n) go to 390 
         if (g  ==  0.0_wp) go to 360 
!                                                                       
         do 320 j = l, n 
!     .......... double division avoids possible underflow ..........   
  320    v(j,i) = (u(i,j) / u(i,l)) / g 
!                                                                       
         do 350 j = l, n 
            s = 0.0_wp 
!                                                                       
            do 340 k = l, n 
  340       s = s + u(i,k) * v(k,j) 
!                                                                       
            do 350 k = l, n 
               v(k,j) = v(k,j) + s * v(k,i) 
  350    continue 
!                                                                       
  360    do 380 j = l, n 
            v(i,j) = 0.0_wp 
            v(j,i) = 0.0_wp 
  380    continue 
!                                                                       
  390    v(i,i) = 1.0_wp 
         g = rv1(i) 
         l = i 
  400 continue 
!     .......... accumulation of left-hand transformations ..........   
  410 if (.not. matu) go to 510 
!     ..........for i=min(m,n) step -1 until 1 do -- ..........         
      mn = n 
      if (m  <  n) mn = m 
!                                                                       
      do 500 ii = 1, mn 
         i = mn + 1 - ii 
         l = i + 1 
         g = w(i) 
         if (i  ==  n) go to 430 
!                                                                       
         do 420 j = l, n 
  420    u(i,j) = 0.0_wp 
!                                                                       
  430    if (g  ==  0.0_wp) go to 475 
         if (i  ==  mn) go to 460 
!                                                                       
         do 450 j = l, n 
            s = 0.0_wp 
!                                                                       
            do 440 k = l, m 
  440       s = s + u(k,i) * u(k,j) 
!     .......... double division avoids possible underflow ..........   
            f = (s / u(i,i)) / g 
!                                                                       
            do 450 k = i, m 
               u(k,j) = u(k,j) + f * u(k,i) 
  450    continue 
!                                                                       
  460    do 470 j = i, m 
  470    u(j,i) = u(j,i) / g 
!                                                                       
         go to 490 
!                                                                       
  475    do 480 j = i, m 
  480    u(j,i) = 0.0_wp 
!                                                                       
  490    u(i,i) = u(i,i) + 1.0_wp 
  500 continue 
!     .......... diagonalization of the bidiagonal form ..........      
  510 tst1 = x 
!     .......... for k=n step -1 until 1 do -- ..........               
      do 700 kk = 1, n 
         k1 = n - kk 
         k = k1 + 1 
         its = 0 
!     .......... test for splitting.                                    
!                for l=k step -1 until 1 do -- ..........               
  520    do 530 ll = 1, k 
            l1 = k - ll 
            l = l1 + 1 
            tst2 = tst1 + abs(rv1(l)) 
            if (tst2  ==  tst1) go to 565 
!     .......... rv1(1) is always zero, so there is no exit             
!                through the bottom of the loop ..........              
            tst2 = tst1 + abs(w(l1)) 
            if (tst2  ==  tst1) go to 540 
  530    continue 
!     .......... cancellation of rv1(l) if l greater than 1 ..........  
  540    c = 0.0_wp 
         s = 1.0_wp 
!                                                                       
         do 560 i = l, k 
            f = s * rv1(i) 
            rv1(i) = c * rv1(i) 
            tst2 = tst1 + abs(f) 
            if (tst2  ==  tst1) go to 565 
            g = w(i) 
            h = pythag(f,g) 
            w(i) = h 
            c = g / h 
            s = -f / h 
            if (.not. matu) go to 560 
!                                                                       
            do 550 j = 1, m 
               y = u(j,l1) 
               z = u(j,i) 
               u(j,l1) = y * c + z * s 
               u(j,i) = -y * s + z * c 
  550       continue 
!                                                                       
  560    continue 
!     .......... test for convergence ..........                        
  565    z = w(k) 
         if (l  ==  k) go to 650 
!     .......... shift from bottom 2 by 2 minor ..........              
         if (its  ==  30) go to 1000 
         its = its + 1 
         x = w(l) 
         y = w(k1) 
         g = rv1(k1) 
         h = rv1(k) 
         f = 0.5_wp* (((g + z) / h) * ((g - z) / y) + y / h - h / y) 
         g = pythag(f,1.0_wp) 
         f = x - (z / x) * z + (h / x) * (y / (f + sign(g,f)) - h) 
!     .......... next qr transformation ..........                      
         c = 1.0_wp 
         s = 1.0_wp 
!                                                                       
         do 600 i1 = l, k1 
            i = i1 + 1 
            g = rv1(i) 
            y = w(i) 
            h = s * g 
            g = c * g 
            z = pythag(f,h) 
            rv1(i1) = z 
            c = f / z 
            s = h / z 
            f = x * c + g * s 
            g = -x * s + g * c 
            h = y * s 
            y = y * c 
            if (.not. matv) go to 575 
!                                                                       
            do 570 j = 1, n 
               x = v(j,i1) 
               z = v(j,i) 
               v(j,i1) = x * c + z * s 
               v(j,i) = -x * s + z * c 
  570       continue 
!                                                                       
  575       z = pythag(f,h) 
            w(i1) = z 
!     .......... rotation can be arbitrary if z is zero ..........      
            if (z  ==  0.0_wp) go to 580 
            c = f / z 
            s = h / z 
  580       f = c * g + s * y 
            x = -s * g + c * y 
            if (.not. matu) go to 600 
!                                                                       
            do 590 j = 1, m 
               y = u(j,i1) 
               z = u(j,i) 
               u(j,i1) = y * c + z * s 
               u(j,i) = -y * s + z * c 
  590       continue 
!                                                                       
  600    continue 
!                                                                       
         rv1(l) = 0.0_wp 
         rv1(k) = f 
         w(k) = x 
         go to 520 
!     .......... convergence ..........                                 
  650    if (z  >=  0.0_wp) go to 700 
!     .......... w(k) is made non-negative ..........                   
         w(k) = -z 
         if (.not. matv) go to 700 
!                                                                       
         do 690 j = 1, n 
  690    v(j,k) = -v(j,k) 
!                                                                       
  700 continue 
!                                                                       
      go to 1001 
!     .......... set error -- no convergence to a                       
!                singular value after 30 iterations ..........          
 1000 ierr = k 
 1001 return 

      end subroutine svd

!=============================================================================80

      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,                          &
     &                  ierr,rv1,rv2,rv3,rv4,rv6)                       
!                                                                       
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group 
      real(wp) d(n),e(n),e2(n),w(m),z(nm,m),                    &
     &       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)                         
      real(wp) u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order
      integer ind(m) 
!                                                                       
!     this subroutine is a translation of the inverse iteration tech-   
!     nique in the algol procedure tristurm by peters and wilkinson.    
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).   
!                                                                       
!     this subroutine finds those eigenvectors of a tridiagonal         
!     symmetric matrix corresponding to specified eigenvalues,          
!     using inverse iteration.                                          
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        d contains the diagonal elements of the input matrix.          
!                                                                       
!        e contains the subdiagonal elements of the input matrix        
!          in its last n-1 positions.  e(1) is arbitrary.               
!                                                                       
!        e2 contains the squares of the corresponding elements of e,    
!          with zeros corresponding to negligible elements of e.        
!          e(i) is considered negligible if it is not larger than       
!          the product of the relative machine precision and the sum    
!          of the magnitudes of d(i) and d(i-1).  e2(1) must contain    
!          0.0_wp if the eigenvalues are in ascending order, or 2.0_wp    
!          if the eigenvalues are in descending order.  if  bisect,     
!          tridib, or  imtqlv  has been used to find the eigenvalues,   
!          their output e2 array is exactly what is expected here.      
!                                                                       
!        m is the number of specified eigenvalues.                      
!                                                                       
!        w contains the m eigenvalues in ascending or descending order. 
!                                                                       
!        ind contains in its first m positions the submatrix indices    
!          associated with the corresponding eigenvalues in w --        
!          1 for eigenvalues belonging to the first submatrix from      
!          the top, 2 for those belonging to the second submatrix, etc. 
!                                                                       
!     on output                                                         
!                                                                       
!        all input arrays are unaltered.                                
!                                                                       
!        z contains the associated set of orthonormal eigenvectors.     
!          any vector which fails to converge is set to zero.           
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          -r         if the eigenvector corresponding to the r-th      
!                     eigenvalue fails to converge in 5 iterations.     
!                                                                       
!        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.      
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      if (m  ==  0) go to 1001 
      tag = 0 
      order = 1.0_wp - e2(1) 
      q = 0 
!     .......... establish and process next submatrix ..........        
  100 p = q + 1 
!                                                                       
      do 120 q = p, n 
         if (q  ==  n) go to 140 
         if (e2(q+1)  ==  0.0_wp) go to 140 
  120 continue 
!     .......... find vectors by inverse iteration ..........           
  140 tag = tag + 1 
      s = 0 
!                                                                       
      do 920 r = 1, m 
         if (ind(r)  /=  tag) go to 920 
         its = 1 
         x1 = w(r) 
         if (s  /=  0) go to 510 
!     .......... check for isolated root ..........                     
         xu = 1.0_wp 
         if (p  /=  q) go to 490 
         rv6(p) = 1.0_wp 
         go to 870 
  490    norm = abs(d(p)) 
         ip = p + 1 
!                                                                       
         do 500 i = ip, q 
  500    norm = max(norm, abs(d(i))+abs(e(i))) 
!     .......... eps2 is the criterion for grouping,                    
!                eps3 replaces zero pivots and equal                    
!                roots are modified by eps3,                            
!                eps4 is taken very small to avoid overflow ..........  
!        eps2 = 1.0d-3 * norm 
         eps2 = 0.001_wp * norm 
         eps3 = epslon(norm) 
         uk = q - p + 1 
         eps4 = uk * eps3 
         uk = eps4 / sqrt(uk) 
         s = p 
  505    group = 0 
         go to 520 
!     .......... look for close or coincident roots ..........          
  510    if (abs(x1-x0)  >=  eps2) go to 505 
         group = group + 1 
         if (order * (x1 - x0)  <=  0.0_wp) x1 = x0 + order * eps3 
!     .......... elimination with interchanges and                      
!                initialization of vector ..........                    
  520    v = 0.0_wp 
!                                                                       
         do 580 i = p, q 
            rv6(i) = uk 
            if (i  ==  p) go to 560 
            if (abs(e(i))  <  abs(u)) go to 540 
!     .......... warning -- a divide check may occur here if            
!                e2 array has not been specified correctly ..........   
            xu = u / e(i) 
            rv4(i) = xu 
            rv1(i-1) = e(i) 
            rv2(i-1) = d(i) - x1 
            rv3(i-1) = 0.0_wp 
            if (i  /=  q) rv3(i-1) = e(i+1) 
            u = v - xu * rv2(i-1) 
            v = -xu * rv3(i-1) 
            go to 580 
  540       xu = e(i) / u 
            rv4(i) = xu 
            rv1(i-1) = u 
            rv2(i-1) = v 
            rv3(i-1) = 0.0_wp 
  560       u = d(i) - x1 - xu * v 
            if (i  /=  q) v = e(i+1) 
  580    continue 
!                                                                       
         if (u  ==  0.0_wp) u = eps3 
         rv1(q) = u 
         rv2(q) = 0.0_wp 
         rv3(q) = 0.0_wp 
!     .......... back substitution                                      
!                for i=q step -1 until p do -- ..........               
  600    do 620 ii = p, q 
            i = p + q - ii 
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i) 
            v = u 
            u = rv6(i) 
  620    continue 
!     .......... orthogonalize with respect to previous                 
!                members of group ..........                            
         if (group  ==  0) go to 700 
         j = r 
!                                                                       
         do 680 jj = 1, group 
  630       j = j - 1 
            if (ind(j)  /=  tag) go to 630 
            xu = 0.0_wp 
!                                                                       
            do 640 i = p, q 
  640       xu = xu + rv6(i) * z(i,j) 
!                                                                       
            do 660 i = p, q 
  660       rv6(i) = rv6(i) - xu * z(i,j) 
!                                                                       
  680    continue 
!                                                                       
  700    norm = 0.0_wp 
!                                                                       
         do 720 i = p, q 
  720    norm = norm + abs(rv6(i)) 
!                                                                       
         if (norm  >=  1.0_wp) go to 840 
!     .......... forward substitution ..........                        
         if (its  ==  5) go to 830 
         if (norm  /=  0.0_wp) go to 740 
         rv6(s) = eps4 
         s = s + 1 
         if (s  >  q) s = p 
         go to 780 
  740    xu = eps4 / norm 
!                                                                       
         do 760 i = p, q 
  760    rv6(i) = rv6(i) * xu 
!     .......... elimination operations on next vector                  
!                iterate ..........                                     
  780    do 820 i = ip, q 
            u = rv6(i) 
!     .......... if rv1(i-1)  ==  e(i), a row interchange               
!                was performed earlier in the                           
!                triangularization process ..........                   
            if (rv1(i-1)  /=  e(i)) go to 800 
            u = rv6(i-1) 
            rv6(i-1) = rv6(i) 
  800       rv6(i) = u - rv4(i) * rv6(i-1) 
  820    continue 
!                                                                       
         its = its + 1 
         go to 600 
!     .......... set error -- non-converged eigenvector ..........      
  830    ierr = -r 
         xu = 0.0_wp 
         go to 870 
!     .......... normalize so that sum of squares is                    
!                1 and expand to full order ..........                  
  840    u = 0.0_wp 
!                                                                       
         do 860 i = p, q 
  860    u = pythag(u,rv6(i)) 
!                                                                       
         xu = 1.0_wp / u 
!                                                                       
  870    do 880 i = 1, n 
  880    z(i,r) = 0.0_wp 
!                                                                       
         do 900 i = p, q 
  900    z(i,r) = rv6(i) * xu 
!                                                                       
         x0 = x1 
  920 continue 
!                                                                       
      if (q  <  n) go to 100 
 1001 return 
      end subroutine tinvit

!=============================================================================80

      subroutine tql1(n,d,e,ierr) 
!                                                                       
      integer i,j,l,m,n,ii,l1,l2,mml,ierr 
      real(wp) d(n),e(n) 
      real(wp) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
!                                                                       
!     this subroutine is a translation of the algol procedure tql1,     
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and     
!     wilkinson.                                                        
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).   
!                                                                       
!     this subroutine finds the eigenvalues of a symmetric              
!     tridiagonal matrix by the ql method.                              
!                                                                       
!     on input                                                          
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        d contains the diagonal elements of the input matrix.          
!                                                                       
!        e contains the subdiagonal elements of the input matrix        
!          in its last n-1 positions.  e(1) is arbitrary.               
!                                                                       
!      on output                                                        
!                                                                       
!        d contains the eigenvalues in ascending order.  if an          
!          error exit is made, the eigenvalues are correct and          
!          ordered for indices 1,2,...ierr-1, but may not be            
!          the smallest eigenvalues.                                    
!                                                                       
!        e has been destroyed.                                          
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the j-th eigenvalue has not been               
!                     determined after 30 iterations.                   
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      if (n  ==  1) go to 1001 
!                                                                       
      do 100 i = 2, n 
  100 e(i-1) = e(i) 
!                                                                       
      f = 0.0_wp 
      tst1 = 0.0_wp 
      e(n) = 0.0_wp 
!                                                                       
      do 290 l = 1, n 
         j = 0 
         h = abs(d(l)) + abs(e(l)) 
         if (tst1  <  h) tst1 = h 
!     .......... look for small sub-diagonal element ..........         
         do 110 m = l, n 
            tst2 = tst1 + abs(e(m)) 
            if (tst2  ==  tst1) go to 120 
!     .......... e(n) is always zero, so there is no exit               
!                through the bottom of the loop ..........              
  110    continue 
!                                                                       
  120    if (m  ==  l) go to 210 
  130    if (j  ==  30) go to 1000 
         j = j + 1 
!     .......... form shift ..........                                  
         l1 = l + 1 
         l2 = l1 + 1 
         g = d(l) 
         p = (d(l1) - g) / (2.0_wp * e(l)) 
         r = pythag(p,1.0_wp) 
         d(l) = e(l) / (p + sign(r,p)) 
         d(l1) = e(l) * (p + sign(r,p)) 
         dl1 = d(l1) 
         h = g - d(l) 
         if (l2  >  n) go to 145 
!                                                                       
         do 140 i = l2, n 
  140    d(i) = d(i) - h 
!                                                                       
  145    f = f + h 
!     .......... ql transformation ..........                           
         p = d(m) 
         c = 1.0_wp 
         c2 = c 
         el1 = e(l1) 
         s = 0.0_wp 
         mml = m - l 
!     .......... for i=m-1 step -1 until l do -- ..........             
         do 200 ii = 1, mml 
            c3 = c2 
            c2 = c 
            s2 = s 
            i = m - ii 
            g = c * e(i) 
            h = c * p 
            r = pythag(p,e(i)) 
            e(i+1) = s * r 
            s = e(i) / r 
            c = p / r 
            p = c * d(i) - s * g 
            d(i+1) = h + s * (c * g + s * d(i)) 
  200    continue 
!                                                                       
         p = -s * s2 * c3 * el1 * e(l) / dl1 
         e(l) = s * p 
         d(l) = c * p 
         tst2 = tst1 + abs(e(l)) 
         if (tst2  >  tst1) go to 130 
  210    p = d(l) + f 
!     .......... order eigenvalues ..........                           
         if (l  ==  1) go to 250 
!     .......... for i=l step -1 until 2 do -- ..........               
         do 230 ii = 2, l 
            i = l + 2 - ii 
            if (p  >=  d(i-1)) go to 270 
            d(i) = d(i-1) 
  230    continue 
!                                                                       
  250    i = 1 
  270    d(i) = p 
  290 continue 
!                                                                       
      go to 1001 
!     .......... set error -- no convergence to an                      
!                eigenvalue after 30 iterations ..........              
 1000 ierr = l 
 1001 return 
      end subroutine tql1

!=============================================================================80

      subroutine tql2(nm,n,d,e,z,ierr) 
!                                                                       
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr 
      real(wp) d(n),e(n),z(nm,n) 
      real(wp) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
!                                                                       
!     this subroutine is a translation of the algol procedure tql2,     
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and     
!     wilkinson.                                                        
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).   
!                                                                       
!     this subroutine finds the eigenvalues and eigenvectors            
!     of a symmetric tridiagonal matrix by the ql method.               
!     the eigenvectors of a full symmetric matrix can also              
!     be found if  tred2  has been used to reduce this                  
!     full matrix to tridiagonal form.                                  
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        d contains the diagonal elements of the input matrix.          
!                                                                       
!        e contains the subdiagonal elements of the input matrix        
!          in its last n-1 positions.  e(1) is arbitrary.               
!                                                                       
!        z contains the transformation matrix produced in the           
!          reduction by  tred2, if performed.  if the eigenvectors      
!          of the tridiagonal matrix are desired, z must contain        
!          the identity matrix.                                         
!                                                                       
!      on output                                                        
!                                                                       
!        d contains the eigenvalues in ascending order.  if an          
!          error exit is made, the eigenvalues are correct but          
!          unordered for indices 1,2,...,ierr-1.                        
!                                                                       
!        e has been destroyed.                                          
!                                                                       
!        z contains orthonormal eigenvectors of the symmetric           
!          tridiagonal (or full) matrix.  if an error exit is made,     
!          z contains the eigenvectors associated with the stored       
!          eigenvalues.                                                 
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the j-th eigenvalue has not been               
!                     determined after 30 iterations.                   
!                                                                       
!     calls pythag for  sqrt(a*a + b*b) .                              
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      ierr = 0 
      if (n  ==  1) go to 1001 
!                                                                       
      do 100 i = 2, n 
  100 e(i-1) = e(i) 
!                                                                       
      f = 0.0_wp 
      tst1 = 0.0_wp 
      e(n) = 0.0_wp 
!                                                                       
      do 240 l = 1, n 
         j = 0 
         h = abs(d(l)) + abs(e(l)) 
         if (tst1  <  h) tst1 = h 
!     .......... look for small sub-diagonal element ..........         
         do 110 m = l, n 
            tst2 = tst1 + abs(e(m)) 
            if (tst2  ==  tst1) go to 120 
!     .......... e(n) is always zero, so there is no exit               
!                through the bottom of the loop ..........              
  110    continue 
!                                                                       
  120    if (m  ==  l) go to 220 
  130    if (j  ==  30) go to 1000 
         j = j + 1 
!     .......... form shift ..........                                  
         l1 = l + 1 
         l2 = l1 + 1 
         g = d(l) 
         p = (d(l1) - g) / (2.0_wp * e(l)) 
         r = pythag(p,1.0_wp) 
         d(l) = e(l) / (p + sign(r,p)) 
         d(l1) = e(l) * (p + sign(r,p)) 
         dl1 = d(l1) 
         h = g - d(l) 
         if (l2  >  n) go to 145 
!                                                                       
         do 140 i = l2, n 
  140    d(i) = d(i) - h 
!                                                                       
  145    f = f + h 
!     .......... ql transformation ..........                           
         p = d(m) 
         c = 1.0_wp 
         c2 = c 
         el1 = e(l1) 
         s = 0.0_wp 
         mml = m - l 
!     .......... for i=m-1 step -1 until l do -- ..........             
         do 200 ii = 1, mml 
            c3 = c2 
            c2 = c 
            s2 = s 
            i = m - ii 
            g = c * e(i) 
            h = c * p 
            r = pythag(p,e(i)) 
            e(i+1) = s * r 
            s = e(i) / r 
            c = p / r 
            p = c * d(i) - s * g 
            d(i+1) = h + s * (c * g + s * d(i)) 
!     .......... form vector ..........                                 
            do 180 k = 1, n 
               h = z(k,i+1) 
               z(k,i+1) = s * z(k,i) + c * h 
               z(k,i) = c * z(k,i) - s * h 
  180       continue 
!                                                                       
  200    continue 
!                                                                       
         p = -s * s2 * c3 * el1 * e(l) / dl1 
         e(l) = s * p 
         d(l) = c * p 
         tst2 = tst1 + abs(e(l)) 
         if (tst2  >  tst1) go to 130 
  220    d(l) = d(l) + f 
  240 continue 
!     .......... order eigenvalues and eigenvectors ..........          
      do 300 ii = 2, n 
         i = ii - 1 
         k = i 
         p = d(i) 
!                                                                       
         do 260 j = ii, n 
            if (d(j)  >=  p) go to 260 
            k = j 
            p = d(j) 
  260    continue 
!                                                                       
         if (k  ==  i) go to 300 
         d(k) = d(i) 
         d(i) = p 
!                                                                       
         do 280 j = 1, n 
            p = z(j,i) 
            z(j,i) = z(j,k) 
            z(j,k) = p 
  280    continue 
!                                                                       
  300 continue 
!                                                                       
      go to 1001 
!     .......... set error -- no convergence to an                      
!                eigenvalue after 30 iterations ..........              
 1000 ierr = l 
 1001 return 
      end subroutine tql2

!=============================================================================80

      subroutine trbak1(nm,n,a,e,m,z) 
!                                                                       
      integer i,j,k,l,m,n,nm 
      real(wp) a(nm,n),e(n),z(nm,m) 
      real(wp) s 
!                                                                       
!     this subroutine is a translation of the algol procedure trbak1,   
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
!                                                                       
!     this subroutine forms the eigenvectors of a real symmetric        
!     matrix by back transforming those of the corresponding            
!     symmetric tridiagonal matrix determined by  tred1.                
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        a contains information about the orthogonal trans-             
!          formations used in the reduction by  tred1                   
!          in its strict lower triangle.                                
!                                                                       
!        e contains the subdiagonal elements of the tridiagonal         
!          matrix in its last n-1 positions.  e(1) is arbitrary.        
!                                                                       
!        m is the number of eigenvectors to be back transformed.        
!                                                                       
!        z contains the eigenvectors to be back transformed             
!          in its first m columns.                                      
!                                                                       
!     on output                                                         
!                                                                       
!        z contains the transformed eigenvectors                        
!          in its first m columns.                                      
!                                                                       
!     note that trbak1 preserves vector euclidean norms.                
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (m  ==  0) go to 200 
      if (n  ==  1) go to 200 
!                                                                       
      do 140 i = 2, n 
         l = i - 1 
         if (e(i)  ==  0.0_wp) go to 140 
!                                                                       
         do 130 j = 1, m 
            s = 0.0_wp 
!                                                                       
            do 110 k = 1, l 
  110       s = s + a(i,k) * z(k,j) 
!     .......... divisor below is negative of h formed in tred1.        
!                double division avoids possible underflow ..........   
            s = (s / a(i,l)) / e(i) 
!                                                                       
            do 120 k = 1, l 
  120       z(k,j) = z(k,j) + s * a(i,k) 
!                                                                       
  130    continue 
!                                                                       
  140 continue 
!                                                                       
  200 return 
      end subroutine trbak1

!=============================================================================80

      subroutine tred1(nm,n,a,d,e,e2) 
!                                                                       
      integer i,j,k,l,n,ii,nm,jp1 
      real(wp) a(nm,n),d(n),e(n),e2(n) 
      real(wp) f,g,h,scale 
!                                                                       
!     this subroutine is a translation of the algol procedure tred1,    
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
!                                                                       
!     this subroutine reduces a real symmetric matrix                   
!     to a symmetric tridiagonal matrix using                           
!     orthogonal similarity transformations.                            
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        a contains the real symmetric input matrix.  only the          
!          lower triangle of the matrix need be supplied.               
!                                                                       
!     on output                                                         
!                                                                       
!        a contains information about the orthogonal trans-             
!          formations used in the reduction in its strict lower         
!          triangle.  the full upper triangle of a is unaltered.        
!                                                                       
!        d contains the diagonal elements of the tridiagonal matrix.    
!                                                                       
!        e contains the subdiagonal elements of the tridiagonal         
!          matrix in its last n-1 positions.  e(1) is set to zero.      
!                                                                       
!        e2 contains the squares of the corresponding elements of e.    
!          e2 may coincide with e if the squares are not needed.        
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      do 100 i = 1, n 
         d(i) = a(n,i) 
         a(n,i) = a(i,i) 
  100 continue 
!     .......... for i=n step -1 until 1 do -- ..........               
      do 300 ii = 1, n 
         i = n + 1 - ii 
         l = i - 1 
         h = 0.0_wp 
         scale = 0.0_wp 
         if (l  <  1) go to 130 
!     .......... scale row (algol tol then not needed) ..........       
         do 120 k = 1, l 
  120    scale = scale + abs(d(k)) 
!                                                                       
         if (scale  /=  0.0_wp) go to 140 
!                                                                       
         do 125 j = 1, l 
            d(j) = a(l,j) 
            a(l,j) = a(i,j) 
            a(i,j) = 0.0_wp 
  125    continue 
!                                                                       
  130    e(i) = 0.0_wp 
         e2(i) = 0.0_wp 
         go to 300 
!                                                                       
  140    do 150 k = 1, l 
            d(k) = d(k) / scale 
            h = h + d(k) * d(k) 
  150    continue 
!                                                                       
         e2(i) = scale * scale * h 
         f = d(l) 
         g = -sign(sqrt(h),f) 
         e(i) = scale * g 
         h = h - f * g 
         d(l) = f - g 
         if (l  ==  1) go to 285 
!     .......... form a*u ..........                                    
         do 170 j = 1, l 
  170    e(j) = 0.0_wp 
!                                                                       
         do 240 j = 1, l 
            f = d(j) 
            g = e(j) + a(j,j) * f 
            jp1 = j + 1 
            if (l  <  jp1) go to 220 
!                                                                       
            do 200 k = jp1, l 
               g = g + a(k,j) * d(k) 
               e(k) = e(k) + a(k,j) * f 
  200       continue 
!                                                                       
  220       e(j) = g 
  240    continue 
!     .......... form p ..........                                      
         f = 0.0_wp 
!                                                                       
         do 245 j = 1, l 
            e(j) = e(j) / h 
            f = f + e(j) * d(j) 
  245    continue 
!                                                                       
         h = f / (h + h) 
!     .......... form q ..........                                      
         do 250 j = 1, l 
  250    e(j) = e(j) - h * d(j) 
!     .......... form reduced a ..........                              
         do 280 j = 1, l 
            f = d(j) 
            g = e(j) 
!                                                                       
            do 260 k = j, l 
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k) 
!                                                                       
  280    continue 
!                                                                       
  285    do 290 j = 1, l 
            f = d(j) 
            d(j) = a(l,j) 
            a(l,j) = a(i,j) 
            a(i,j) = f * scale 
  290    continue 
!                                                                       
  300 continue 
!                                                                       
      return 
      end subroutine tred1

!=============================================================================80

      subroutine tred2(nm,n,a,d,e,z) 
!                                                                       
      integer i,j,k,l,n,ii,nm,jp1 
      real(wp) a(nm,n),d(n),e(n),z(nm,n) 
      real(wp) f,g,h,hh,scale 
!                                                                       
!     this subroutine is a translation of the algol procedure tred2,    
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
!                                                                       
!     this subroutine reduces a real symmetric matrix to a              
!     symmetric tridiagonal matrix using and accumulating               
!     orthogonal similarity transformations.                            
!                                                                       
!     on input                                                          
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement.                                         
!                                                                       
!        n is the order of the matrix.                                  
!                                                                       
!        a contains the real symmetric input matrix.  only the          
!          lower triangle of the matrix need be supplied.               
!                                                                       
!     on output                                                         
!                                                                       
!        d contains the diagonal elements of the tridiagonal matrix.    
!                                                                       
!        e contains the subdiagonal elements of the tridiagonal         
!          matrix in its last n-1 positions.  e(1) is set to zero.      
!                                                                       
!        z contains the orthogonal transformation matrix                
!          produced in the reduction.                                   
!                                                                       
!        a and z may coincide.  if distinct, a is unaltered.            
!                                                                       
!     questions and comments should be directed to burton s. garbow,    
!     mathematics and computer science div, argonne national laboratory 
!                                                                       
!     this version dated august 1983.                                   
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      do 100 i = 1, n 
!                                                                       
         do 80 j = i, n 
   80    z(j,i) = a(j,i) 
!                                                                       
         d(i) = a(n,i) 
  100 continue 
!                                                                       
      if (n  ==  1) go to 510 
!     .......... for i=n step -1 until 2 do -- ..........               
      do 300 ii = 2, n 
         i = n + 2 - ii 
         l = i - 1 
         h = 0.0_wp 
         scale = 0.0_wp 
         if (l  <  2) go to 130 
!     .......... scale row (algol tol then not needed) ..........       
         do 120 k = 1, l 
  120    scale = scale + abs(d(k)) 
!                                                                       
         if (scale  /=  0.0_wp) go to 140 
  130    e(i) = d(l) 
!                                                                       
         do 135 j = 1, l 
            d(j) = z(l,j) 
            z(i,j) = 0.0_wp 
            z(j,i) = 0.0_wp 
  135    continue 
!                                                                       
         go to 290 
!                                                                       
  140    do 150 k = 1, l 
            d(k) = d(k) / scale 
            h = h + d(k) * d(k) 
  150    continue 
!                                                                       
         f = d(l) 
         g = -sign(sqrt(h),f) 
         e(i) = scale * g 
         h = h - f * g 
         d(l) = f - g 
!     .......... form a*u ..........                                    
         do 170 j = 1, l 
  170    e(j) = 0.0_wp 
!                                                                       
         do 240 j = 1, l 
            f = d(j) 
            z(j,i) = f 
            g = e(j) + z(j,j) * f 
            jp1 = j + 1 
            if (l  <  jp1) go to 220 
!                                                                       
            do 200 k = jp1, l 
               g = g + z(k,j) * d(k) 
               e(k) = e(k) + z(k,j) * f 
  200       continue 
!                                                                       
  220       e(j) = g 
  240    continue 
!     .......... form p ..........                                      
         f = 0.0_wp 
!                                                                       
         do 245 j = 1, l 
            e(j) = e(j) / h 
            f = f + e(j) * d(j) 
  245    continue 
!                                                                       
         hh = f / (h + h) 
!     .......... form q ..........                                      
         do 250 j = 1, l 
  250    e(j) = e(j) - hh * d(j) 
!     .......... form reduced a ..........                              
         do 280 j = 1, l 
            f = d(j) 
            g = e(j) 
!                                                                       
            do 260 k = j, l 
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k) 
!                                                                       
            d(j) = z(l,j) 
            z(i,j) = 0.0_wp 
  280    continue 
!                                                                       
  290    d(i) = h 
  300 continue 
!     .......... accumulation of transformation matrices ..........     
      do 500 i = 2, n 
         l = i - 1 
         z(n,l) = z(l,l) 
         z(l,l) = 1.0_wp 
         h = d(i) 
         if (h  ==  0.0_wp) go to 380 
!                                                                       
         do 330 k = 1, l 
  330    d(k) = z(k,i) / h 
!                                                                       
         do 360 j = 1, l 
            g = 0.0_wp 
!                                                                       
            do 340 k = 1, l 
  340       g = g + z(k,i) * z(k,j) 
!                                                                       
            do 360 k = 1, l 
               z(k,j) = z(k,j) - g * d(k) 
  360    continue 
!                                                                       
  380    do 400 k = 1, l 
  400    z(k,i) = 0.0_wp 
!                                                                       
  500 continue 
!                                                                       
  510 do 520 i = 1, n 
         d(i) = z(n,i) 
         z(n,i) = 0.0_wp 
  520 continue 
!                                                                       
      z(n,n) = 1.0_wp 
      e(1) = 0.0_wp 
      return 
      end subroutine tred2

!=============================================================================80
      subroutine Eigen_Sym(n,AA,RR,EE)

      implicit none

      integer,                     intent(in   ) :: n
      real(wp), dimension(n,n), intent(in   ) :: AA

      real(wp), dimension(n  ), intent(  out) :: EE
      real(wp), dimension(n,n), intent(  out) :: RR

!     input
!     n : dimension of matrix
!     A : AA(n,n) ->  A = A^T

!     output
!     R   : Matrix of eigenvectors satisfying  A = R lam R^T ; 
!           with     R R^T = R R^T = I
!     Eig : Vector of eigenvalues

!     Determine the eigenvalues and eigenvectors for a symmetric matrix
!     we do not assume the matrix A is SPD, only symmetric
!     
      real(wp), parameter                     :: eps = 1.0e-15_wp
      real(wp) :: root
      real(wp) :: t11,t12,t21,t22
      real(wp) :: a11,a12,a13,a22,a23,a33
      real(wp) :: b,c,d, p,q
      real(wp) :: mag1,mag2
      real(wp) :: pi, disc, phi
      real(wp) :: r1,r2,r3,t1, err, mag, al
      complex(wp) :: z3p,z3m
      integer        :: j

      continue

      pi = acos(-1.0_wp)
      if(maxval(abs(AA - Transpose(AA))) >= eps) then
        write(*,*)'matrix is not symmetric;  stopping'
      endif

      if(n == 1) then

        write(*,*)'matrix dimension == 1;  stopping'
        stop

      elseif(n == 2) then

        a11 = AA(1,1) ;  a12 = AA(1,2)
                         a22 = AA(2,2) 

        root   = sqrt(4.0_wp*a12*a12 +(a11 - a22)*(a11 - a22))
        EE(1) = 0.5_wp * ((a11 + a22) - root )
        EE(2) = 0.5_wp * ((a11 + a22) + root )

        t11 = -( -a11 + a22 + root )/(2.0_wp*a12)
        t12 = 1.0_wp
        mag1 = 1.0_wp / sqrt(t11*t11+t12*t12)
        
        t21 = +( +a11 - a22 + root )/(2.0_wp*a12)
        t22 = 1.0_wp
        mag2 = 1.0_wp / sqrt(t21*t21+t22*t22)

        RR(1,1) = t11 * mag1
        RR(1,2) = t12 * mag1
        RR(2,1) = t21 * mag2
        RR(2,2) = t12 * mag2

      elseif(n == 3) then

         a11 = AA(1,1) ; a12 = AA(1,2) ; a13 = AA(1,3)
                         a22 = AA(2,2) ; a23 = AA(2,3)
                                         a33 = AA(3,3)

         b = (- a11 - a22 - a33)
         c = (- a12*a12 - a13*a13 + a11*a22                 &
           &  - a23*a23 + a11*a33 + a22*a33)
         d = (+ a13*a13*a22 - a12*a13*a23 - a12*a13*a23     &
           &  + a11*a23*a23 + a12*a12*a33 - a11*a22*a33)

         p  = c - b*b / 3.0_wp
         q  = d + (2.0_wp*b*b*b - 9.0_wp*b*c)/(27.0_wp)

         phi = - acos(1.5_wp*q/p*sqrt(-3.0_wp/p))
         al = 2.0_wp * sqrt(-p/3.0_wp) 
         EE(1) =  al * cos((phi - 2.0_wp*0.0_wp*pi)/3.0_wp) - b/3.0_wp
         EE(2) =  al * cos((phi - 2.0_wp*1.0_wp*pi)/3.0_wp) - b/3.0_wp
         EE(3) =  al * cos((phi - 2.0_wp*2.0_wp*pi)/3.0_wp) - b/3.0_wp
!        err =                                                     &
!            & + EE(1)*EE(1)*EE(1) + b*EE(1)*EE(1) + c*EE(1) + d   &
!            & + EE(2)*EE(2)*EE(2) + b*EE(2)*EE(2) + c*EE(2) + d   &
!            & + EE(3)*EE(3)*EE(3) + b*EE(3)*EE(3) + c*EE(3) + d
!        if(sqrt(err*err) >= 1.0e-11) then
!          write(*,*)'cubic polynomial error:  stopping'
!          stop
!        endif
         do j = 1,3
           RR(1,j) = 1.0_wp  
           RR(2,j) = - (a13*a23 + a12*(-a33 + EE(j))) / (a23*a23 + (a22 - EE(j))*(-a33 + EE(j)))
           RR(3,j) = - (a12*a23 + a13*(-a22 + EE(j))) / (a23*a23 + (a22 - EE(j))*(-a33 + EE(j)))
           mag = sqrt(dot_product(RR(:,j),RR(:,j)))
           RR(:,j) = RR(:,j)/mag
         enddo

      elseif(n == 4) then

      elseif(n == 5) then

      else

      endif

      end subroutine Eigen_Sym

      end module eispack_module
