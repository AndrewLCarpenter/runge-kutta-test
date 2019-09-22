pts = 65 ;
mpt = (pts-1)/2 ;
rhs  = zeros(1,pts) ;
vecS = zeros(pts,pts) ;
vecD1= zeros(pts,pts) ;
vecD2= zeros(pts,pts) ;
DFT = zeros(pts,pts) ;
vecS(:,1) = 1 ;
vecD1(:,1) = 0 ;
vecD2(:,1) = 0 ;
DFT(:,1) = 1*sqrt(1/pts) ;

for j = 1:mpt
  for i = 1:pts
    x = 2*pi*(i-1)/(pts+0) ;
   vecS(i,2*j-0) =     sin(j*x) ;
   vecS(i,2*j+1) =     cos(j*x) ;
    DFT(i,2*j-0)  =    sin(j*x)*sqrt(2/pts) ;
    DFT(i,2*j+1)  =    cos(j*x)*sqrt(2/pts) ;
  vecD1(i,2*j-0) =  +j*cos(j*x) ;
  vecD2(i,2*j-0) =-j*j*sin(j*x) ;
  vecD1(i,2*j+1) =  -j*sin(j*x) ;
  vecD2(i,2*j+1) =-j*j*cos(j*x) ;
  end
end
Dmat1 = vecD1*inv(vecS) ;
Dmat2 = Dmat1*Dmat1 ;
Dmat4 = Dmat2*Dmat2 ;
[uu,dd] = schur(Dmat1) ;

% FD operators
ncol2 = [0,-1,+0,0,0,0,0,-0,+1]/2;
ncol4 = [0,-8,+1,0,0,0,0,-1,+8]/12;
ncol6 = [0,-45,+9,-1,0,0,+1,-9,+45]/60;
ncol8 = [0,-672,+168,-32,+3,-3,+32,-168,+672]/840;


NN = pts ;
ncol = zeros(NN,1) ;
dx = 2*pi/NN ;
% 8th-order
ncol(   5) = +   3/840 ;
ncol(   4) = -  32/840 ;
ncol(   3) = + 168/840 ;
ncol(   2) = - 672/840 ;
ncol(NN-0) = + 672/840 ;
ncol(NN-1) = - 168/840 ;
ncol(NN-2) = +  32/840 ;
ncol(NN-3) = -   3/840 ;
% 2nd-order
%ncol = zeros(NN,1) ;
%ncol(   2) = - 1/2 ;
%ncol(NN-0) = + 1/2 ;
nrow = -ncol ;
D2nd = toeplitz(ncol,nrow)/dx ;
Deig = sort(eig(D2nd))   

%Len = 7 ;
%  q = 2*pi/Len ;
%vecT = zeros(pts,pts) ;

%for j = 1:mpt
%  for i = 1:pts
%    xx=  (i-1) *Len / (pts+0)  ;
%    k = q*j ;
%   vecT(i,2*j-0) =     sin(k*xx) ;
%   vecT(i,2*j+1) =     cos(k*xx) ;
%  vecD1(i,2*j-0) =  +k*cos(k*xx) ;
%  vecD2(i,2*j-0) =-k*k*sin(k*xx) ;
%  vecD1(i,2*j+1) =  -k*sin(k*xx) ;
%  vecD2(i,2*j+1) =-k*k*cos(k*xx) ;
%  end
%end
%Dmat1T = vecD1*inv(vecT) ;
%Dmat2T = Dmat1T*Dmat1T ;
%
%eig(Dmat2*q^2 + Dmat4*q^4)


%ncol = zeros(1,pts) ;
%ncol(1)   = -2 ;
%ncol(2)   = +1 ;
%ncol(pts) = +1 ;
%dx = 2*pi/pts ;
%Mmat = toeplitz(ncol) /dx/dx ;
%eig(Mmat + Mmat*Mmat) ;
%L = 8;
%tmp = Dmat2/L^2 + Dmat4/L^4 ;
%eig(tmp) 
%[uu,lam] = eig(tmp) ;

%[U,D,V] = svd(vecs) ;
%Dinv = inv(D) ;
%for i = mpt:pts
%  Dinv(i,i) = 0 ;
%end

%ans = rhs*V*Dinv*U' ;

format long e

Dmat1(:,1)
Dmat2(:,1)
Dmat4(:,1)

