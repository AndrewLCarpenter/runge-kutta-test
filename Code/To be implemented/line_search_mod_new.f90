      module line_search_mod
      
      use precision_vars

      !implicit none
      
      private
      public :: newt
      
      contains

!==============================================================================

      FUNCTION fmin(x)
      INTEGER n,NP
      REAL fmin,x(*),fvec
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
      !  USES funcv
!  Returns f = 1/2 F F at x. subroutine funcv(n,x,f) is  a   xed-name,  user-supplied
!routine that returns the vector of functions at x.  The common block newtv
!communicates the  function  values  back  to newt.
      INTEGER i
      REAL sum
      call funcv(n,x,fvec)
      sum=0.
      do  i=1,n
        sum=sum+fvec(i)**2
      enddo 
        fmin=0.5*sum
      return
      END function fmin     
!==============================================================================
      SUBROUTINE newt(x,n,check)
      INTEGER n,nn,NP,MAXITS
      LOGICAL check
      REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX
      PARAMETER (NP=40,MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7, STPMX=100.)
      COMMON /newtv/ fvec(NP),nn !Communicates with fmin.
      SAVE /newtv/

! USES fdjac,fmin,lnsrch,lubksb,ludcmp
!Given an initial guess x(1:n) for a root in n dimensions, find  the  root  by  a  globally
!convergent  Newton's  method.  The  vector  of  functions  to  be  zeroed,  called
!fvec(1:n) in  the  routine  below,  is  returned  by  a  user-supplied  subroutine  that
!must be  called funcv and  have  the  declaration subroutine funcv(n,x,fvec).  The  output  quantity
!check is false on a normal return and true if the routine has converged to a local minimum of the
!function fmin defined below.  In this  case try restarting from a different  initial  guess.
!Parameters: NP is  the maximum  expected value  of n; MAXITS is the maximum  number of
!iterations; TOLF sets the convergence criterion on function values; TOLMIN
!sets the criterion for  deciding  whether  spurious  convergence  to  a  minimum  of
!fmin has  occurred; TOLX is 380 Chapter 9. Root Finding and Nonlinear Sets of Equations
!Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
!Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by Numerical Recipes Software. 
!Permission is granted for internet users to make one paper copy for their own personal use. Further reproduction, or any copying of machine-
!readable files (including this one) to any server computer, is strictly prohibited. To order Numerical Recipes books, diskettes, or CDROMs
!visit website http://www.nr.com or call 1-800-872-7423 (North America only), or send email to trade@cup.cam.ac.uk (outside North America).
!the convergence criterion on x; STPMX is  the  scaled maximum  step length  allowed  in  line searches.

      INTEGER i,its,j,indx(NP)
      REAL d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP),xold(NP),fmin
      real(wp) :: wrk
      EXTERNAL fmin

      test = maxval(abs(fvec(:)))

      if(test .lt. 0.01_wp*TOLF)then ; check=.false. ; return ; endif ;

        wrk  =  sqrt(dot_product(x(:),x(:)))
      stpmax = STPMX*max(wrk ,float(n))

      do  its=1,MAXITS !Start of  iteration loop.

        call fdjac(n,x,fvec,NP,fjac)

        !  direction of steepest descent  grad(f)
        do  i=1,n !Compute r f for  the line  search.
          g(i) = matmul(Transpose(fjac),fvec)
        enddo 

        xold(:) =      x(:)
           p(:) = - fvec(:)   !  Right-hand side for linear equations.
           fold =   f         !  And f.

        call ludcmp(fjac,n,NP,indx,d) !Solve  linear equations  by LU decomposition.
        call lubksb(fjac,n,NP,indx,p)

        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)

        !  lnsrch returns new x and f.  It  also calculates fvec at the new
        !  x when it  calls fmin.

        test = maxval(abs(fvec(:)))    !Test for convergence on function values.
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then !Check for gradient of f zero,  i.e.,  spurious  convergence.
          test =0.0_wp
           den =max(f,0.5_wp*n)
          do  i=1,n
            temp=abs(g(i))*max(abs(x(i)),1.0_wp) / den
            if(temp.gt.test)test=temp
          enddo 
          if(test.lt.TOLMIN)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=0.0_wp !Test for  convergence on dx.
        do  i=1,n
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.0_wp)
          if(temp.gt.test)test=temp
        enddo 
        if(test.lt.TOLX)return

      enddo 
      WRITE(*,*) 'MAXITS exceeded in newt'
      END subroutine newt

!==============================================================================
!     FUNCTION fmin(x)
!     INTEGER n,NP
!     REAL fmin,x(*),fvec
!     PARAMETER (NP=40)
!     COMMON /newtv/ fvec(NP),n
!     SAVE /newtv/
!     !USES funcv
!     !Returns f = 1 F · F at x. subroutine funcv(n,x,f) is a fixed-name, user-supplied 2
!     !routine that returns the vector of functions at x. The common block newtv communicates
!     !the function values back to newt.
!     INTEGER i
!     REAL sum
!     call funcv(n,x,fvec)
!     sum=0.
!     do  i=1,n
!     sum=sum+fvec(i)**2
!     enddo 
!       fmin=0.5*sum
!     return
!     END function fmin

!==============================================================================
      SUBROUTINE fdjac(n,x,fvec,np,df)
      INTEGER n,np,NMAX
      REAL df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=40,EPS=1.e-4)
      !  USES funcv
!Computes  forward-di erence  approximation  to  Jacobian.  On  input,x(1:n)
!is  the  point at  which  the  Jacobian  is  to  be  evaluated, fvec(1:n)
!is  the  vector  of  function  values  at the  point,  and np
!is  the physical  dimension of  the Jacobian array df(1:n,1:n) which  is output.
!subroutine funcv(n,x,f) is  a   xed-name,  user-supplied  routine  that  returns
!the  vector  of  functions  at x. Parameters: NMAX is  the  maximum  value  of
!n; EPS is  the approximate  square root  of the machine  precision.
     INTEGER i,j
     REAL h,temp,f(NMAX)
      do  j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.)h=EPS
        x(j)=temp+h !Trick  to reduce finite precision error.
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do  i=1,n !Forward difference formula.
          df(i,j)=(f(i)-fvec(i))/h
        enddo 
      enddo 
      return
      END subroutine fdjac

!==============================================================================

      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)

      integer,                intent(in  )  :: n

      REAL(wp),               intent(in   ) :: fold,stpmax
      REAL(wp), dimension(n), intent(in   ) ::  g,p,xold
      REAL(wp), dimension(n), intent(  out) ::  x

      LOGICAL,               intent(inout)  :: check

      integer,  parameter                   :: maxiter = 10
      REAL(wp), parameter                   :: alf=1.0e-4_wp, tolx=1.e-7_wp
      REAL(wp), parameter                   :: tolerance = 1.0e-15_wp

!Given  an n-dimensional  point xold(1:n),  the  value  of  the  function  and  gradient  there,
!fold and g(1:n), and a direction p(1:n),  finds  a new point x(1:n) along the direction
!p from xold where the function funchas decreased \suciently." The new function value
!is returned in f. stpmax is an input quantity that limits the length of the steps so that you
!do not  try to evaluate the function  in  regions where it is  unde ned  or subject to overflow.
!p is  usually  the  Newton  direction.  The  output  quantity check is  false  on  a  normal  exit.
!It  is  true  when x is  too  close  to xold.  In  a  minimization  algorithm,  this  usually  signals
!convergence and can be  ignored.  However, in  a zero- nding  algorithm the  calling  program
!should  check  whether  the  convergence  is  spurious. Parameters: alf ensures  sufficient  decrease in function value;
!tolx is  the  convergence criterion  on x.

      integer               :: i
      real(wp)              :: a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,temp,test,tmplam

      check=.false.

      wrk = sqrt(dot_product(p(:),p(:)))     !  Scale step if it is too big.
      if(wrk.gt.stpmax) p(:) = p(:) * stpmax / wrk    

      slope = dot_product(g(:),p(:))
      if(slope.ge.0.0_wp) write(*,*) 'roundoff problem in lnsrch'
      test=0.0_wp !Compute min.
      do i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.0_wp)
        if(temp.gt.test)test=temp
      enddo 
      alamin= TOLX/test
      alam  = 1.0_wp   ! Always try full Newton step first.

      do k = 1,maxiter ! Start of iteration loop.

        x(:) = xold(:) + alam*p(:)

        f=func(x)

        if(alam.lt.alamin)then !Convergence  on x.  For zero finding, the calling program should verify the convergence.
          x(:)=xold(:)
          check=.true.
          exit
        else if(f.le.fold+ALF*alam*slope)then !Sufficient function decrease.
          exit
        else !Backtrack.
          if(alam.eq.1.0_wp)then !First  time.
            tmplam=-slope/(2.0_wp*(f-fold-slope))
          else !subsequent backtracks.
            rhs1=f -fold-alam *slope
            rhs2=f2-fold-alam2*slope
            a=(       rhs1/alam**2-     rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(abs(a) .le. tolerance)then
              tmplam=-slope/(2.0_wp*b)
            else
              disc=b*b-3.0_wp*a*slope
              if(disc.lt.0.0_wp)then
                tmplam=0.5_wp*alam
              else if(b.le.0.0_wp)then
                tmplam=(-b+sqrt(disc))/(3.0_wp*a)
              else
                tmplam=-slope/(b+sqrt(disc))
              endif
            endif
            if(tmplam .gt. 0.5_wp*alam)tmplam=0.5_wp*alam
          endif
        endif
        alam2=alam
        f2=f
        alam=max(tmplam,0.1_wp*alam)
      enddo
      return

      END subroutine lnsrch

!==============================================================================
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
!Solves the  set of n linear  equations A X = B.Here a is input,  not  as the  matrix
!A but rather  as its LU decomposition, determined by  the routine ludcmp. indx
!is  input  as  the permutation vector returned by ludcmp . b(1:n)is input  as the right-hand side vector
!B, and returns with the solution vector X. a, n , np,and indx are not modi ed by this routine
!and  can  be left  in  place  for successive calls  with  di erent right-hand  sides
!b.  This  routine takes into account the possibility that b will begin with many zero elements, so it is ecient
!for  use  in  matrix  inversion.
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      !When ii is  set  to  a  positive  value,  it  will  become  the  in-
      !dex  of  the   rst  nonvanishing  element  of b.Wenowdo the forward 
      !substitution, equation (2.3.6).  The only new wrinkle  is to unscramble 
      ! the permutation as we go.
      do  i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
        do  j=ii,i-1
          sum=sum-a(i,j)*b(j)
        enddo 
        else if (sum.ne.0.) then
          ii=i
                !A nonzero element was encountered, so from now on we will
                !have  to do  the sums  in  the  loop above.
        endif
        b(i)=sum
      enddo 
      do  i=n,1,-1 !Now we do the backsubstitution, equation  (2.3.7).
        sum=b(i)
        do  j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo 
        b(i)=sum/a(i,i) !Store a  component of the solution vector X.
      enddo 
      return !All  done!
      END subroutine lubksb
!==============================================================================
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20) !Largest expected n, and a small number.
!Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
!the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
!arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
!row permutation effected by the partial pivoting; d is output as ±1 depending on whether
!the number of row interchanges was even or odd, respectively. This routine is used in
!combination with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
!vv stores the implicit scaling of each row. 
      d=1. !No row interchanges yet.
      do  i=1,n !Loop over rows to get the implicit scaling informa-tion.
        aamax=0.
        do  j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo 
         if (aamax.eq.0.) write(*,*),'singular matrix in ludcmp'!No nonzero largest element.
         vv(i)=1./aamax !Save the scaling.
      enddo 
      do  j=1,n !This is the loop over columns of Crout’s method.
        do  i=1,j-1 !This is equation (2.3.12) except for i = j.
          sum=a(i,j)
          do  k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo 
          a(i,j)=sum
        enddo 
        aamax=0. !Initialize for the search for largest pivot element.
        do  i=j,n !This is i = j of equation (2.3.12) and i = j + 1 . . . N
          sum=a(i,j) !of equation (2.3.13).
          do  k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo 
          a(i,j)=sum
          dum=vv(i)*abs(sum) !Figure of merit for the pivot.
          if (dum.ge.aamax) then !Is it better than the best so far?
            imax=i
            aamax=dum
          endif
        enddo 
        if (j.ne.imax)then !Do we need to interchange rows?
          do  k=1,n !Yes, do so...
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo 
          d=-d !...and change the parity of d.
          vv(imax)=vv(j) !Also interchange the scale factor.
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do  i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo 
        endif
      enddo 
      return
      END subroutine ludcmp
!==============================================================================
      SUBROUTINE funcv(n2,v,f)
      INTEGER n2,nvar,kmax,kount,KMAXX,NMAX
      REAL f(n2),v(n2),x1,x2,dxsav,xp,yp,EPS
      PARAMETER (NMAX=50,KMAXX=200,EPS=1.e-6) !At most NMAX coupled ODEs.
      COMMON /caller/ x1,x2,nvar
      COMMON /path/ kmax,kount,dxsav,xp(KMAXX),yp(NMAX,KMAXX)
!USES derivs,load,odeint,rkqs,score
!Routine for use with newt to solve a two point boundary value problem for nvar coupled
!ODEs by shooting from x1 to x2. Initial values for the nvar ODEs at x1 are generated
!from the n2 input coefficients v(1:n2), using the user-supplied routine load. The routine
!integrates the ODEs to x2 using the Runge-Kutta method with tolerance EPS, initial stepsize
!h1, and minimum stepsize hmin. At x2 it calls the user-supplied subroutine score to
!evaluate the n2 functions f(1:n2) that ought to be zero to satisfy the boundary conditions
!at x2. The functions f are returned on output. newt uses a globally convergent Newton’s
!method to adjust the values of v until the functions f are zero. The user-supplied subroutine
!derivs(x,y,dydx) supplies derivative information to the ODE integrator (see Chapter
!16). The common block caller receives its values from the main program so that funcv
!can have the syntax required by newt. The common block path is included for compatibility
!with odeint.
      INTEGER nbad,nok
      REAL h1,hmin,y(NMAX)
      EXTERNAL derivs,rkqs
      kmax=0
      h1=(x2-x1)/100.
      hmin=0.
      call load(x1,v,y)
      call odeint(y,nvar,x1,x2,EPS,h1,hmin,nok,nbad,derivs,rkqs)
      call score(x2,y,f)
      return
      END subroutine funcv
!==============================================================================
      SUBROUTINE derivs(x,y,dydx)
      INTEGER m,n
      REAL c2,dx,gamma,x,dydx(3),y(3)
      COMMON /sphcom/ c2,gamma,dx,m,n
      !Evaluates derivatives for odeint.
      dydx(1)=y(2)
      dydx(2)=(2.0*x*(m+1.0)*y(2)-(y(3)-c2*x*x)*y(1))/(1.0-x*x)
      dydx(3)=0.0
      return
      END subroutine derivs
!==============================================================================
      SUBROUTINE load(x1,v,y)
      INTEGER m,n
      REAL c2,dx,gamma,x1,y1,v(1),y(3)
      COMMON /sphcom/ c2,gamma,dx,m,n
      !Supplies starting values for integration at x = −1 + dx.
      y(3)=v(1)
      if(mod(n-m,2).eq.0)then
        y1=gamma
      else
        y1=-gamma
      endif
      y(2)=-(y(3)-c2)*y1/(2*(m+1))
      y(1)=y1+y(2)*dx
      return
      END subroutine load
!==============================================================================
      SUBROUTINE score(x2,y,f)
      INTEGER m,n
      REAL c2,dx,gamma,x2,f(1),y(3)
      COMMON /sphcom/ c2,gamma,dx,m,n
      !Tests whether boundary condition at x = 0 is satisfied.
      if (mod(n-m,2).eq.0) then
        f(1)=y(2)
      else
        f(1)=y(1)
      endif
      return
      END subroutine score
!==============================================================================
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50) !Maximum number of equations.
!USES derivs,rkck
!Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
!and adjust stepsize. Input are the dependent variable vector y(1:n) and its derivative
!dydx(1:n) at the starting value of the independent variable x. Also input are the stepsize
!to be attempted htry, the required accuracy eps, and the vector yscal(1:n) against
!which the error is scaled. On output, y and x are replaced by their new values, hdid is the
!stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs
!is the user-supplied subroutine that computes the right-hand side derivatives.
      INTEGER i
      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      !The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.
      h=htry !Set stepsize to the initial trial value.
 1    call rkck(y,dydx,n,x,h,ytemp,yerr,derivs) !Take a step.
      errmax=0. !Evaluate accuracy.
      do  i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
      enddo 
      errmax=errmax/eps
      !Scale relative to required tolerance.
      if(errmax.gt.1.)then !Truncation error too large, reduce stepsize.
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)!  No more than a factor of 10.
        xnew=x+h
        if(xnew.eq.x) write(*,*) 'stepsize underflow in rkqs'
        goto 1 !For another try.
      else ! Step succeeded. Compute size of next step.
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else! No more than a factor of 5 increase.
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do  i=1,n
          y(i)=ytemp(i)
        enddo 
        return
      endif
      END subroutine rkqs
!==============================================================================
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
!Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar)
!from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
!h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can
!be zero). On output nok and nbad are the number of good and bad (but retried and
!fixed) steps taken, and ystart is replaced by values at the end of the integration interval.
!derivs is the user-supplied subroutine for calculating the right-hand side derivative, while
!rkqs is the name of the stepper routine to be used. /path/ contains its own information
!about how often an intermediate value is to be stored.
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp !User storage for intermediate results. Preset dxsav and kmax.
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do  i=1,nvar
        y(i)=ystart(i)
      enddo 
        if (kmax.gt.0) xsav=x-2.*dxsav !Assures storage of first step.
        do  nstp=1,MAXSTP !Take at most MAXSTP steps.
          call derivs(x,y,dydx)
          do  i=1,nvar
!Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
            yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
          enddo 
          if(kmax.gt.0)then
            if(abs(x-xsav).gt.abs(dxsav)) then !Store intermediate results.
              if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do  i=1,nvar
                yp(i,kount)=y(i)
              enddo 
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x !If stepsize can overshoot, decrease.
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then !Are we done?
          do  i=1,nvar
            ystart(i)=y(i)
          enddo 
          if(kmax.ne.0)then
            kount=kount+1 !Save final step.
            xp(kount)=x
            do  i=1,nvar
              yp(i,kount)=y(i)
            enddo 
          endif
          return !Normal exit.
        endif
        if(abs(hnext).lt.hmin) write(*,*) 'stepsize smaller than minimum in odeint'
        h=hnext
      enddo  
      write(*,*) 'too many steps in odeint'
      return
      END subroutine odeint
!==============================================================================
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50) !Set to the maximum number of functions.
!USES derivs
!Given values for n variables y and their derivatives dydx known at x, use the fifth-order
!Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
!the incremented variables as yout. Also return an estimate of the local truncation er-
!ror in yout using the embedded fourth-order method. The user supplies the subroutine
!derivs(x,y,dydx), which returns derivatives dydx at x.
      INTEGER i
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX)
      REAL ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51
      REAL B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3
      REAL DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.)
      PARAMETER (B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5)
      PARAMETER (B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.)
      PARAMETER (B63=575./13824.,B64=44275./110592.,B65=253./4096.)
      PARAMETER (C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.)
      PARAMETER (DC1=C1-2825./27648.,DC3=C3-18575./48384.)
      PARAMETER (DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)
      do  i=1,n !First step.
        ytemp(i)=y(i)+B21*h*dydx(i)
      enddo 
      call derivs(x+A2*h,ytemp,ak2) !Second step.
      do  i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      enddo 
      call derivs(x+A3*h,ytemp,ak3) !Third step.
      do  i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      enddo 
      call derivs(x+A4*h,ytemp,ak4) !Fourth step.
      do  i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      enddo 
      call derivs(x+A5*h,ytemp,ak5) !Fifth step.
      do  i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
      enddo 
      call derivs(x+A6*h,ytemp,ak6) !Sixth step.
      do  i=1,n !Accumulate increments with proper weights.
      yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
        enddo 
      do  i=1,n !Estimate error as difference between fourth and fifth order methods.
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
      enddo 
      return
      END subroutine rkck

!==============================================================================
      end module line_search_mod
