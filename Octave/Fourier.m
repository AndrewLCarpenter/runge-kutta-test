pts = 65 ;
mpt = (pts-1)/2 ;
rhs  = zeros(1,pts) ;
vecS = zeros(pts,pts) ;
vecD = zeros(pts,pts) ;
vecS(:,1) = 1 ;
vecD(:,1) = 0 ;

for j = 1:mpt
  for i = 1:pts
    x = 2*pi*(i-1)/(pts+0) ;
  vecS(i,2*j-0) =   sin(j*x) ;
  vecD(i,2*j-0) =+j*cos(j*x) ;
  vecS(i,2*j+1) =   cos(j*x) ;
  vecD(i,2*j+1) =-j*sin(j*x) ;
  end
end

Fourier = vecD*inv(vecS) ;

%[U,D,V] = svd(vecs) ;
%Dinv = inv(D) ;
%for i = mpt:pts
%  Dinv(i,i) = 0 ;
%end

%ans = rhs*V*Dinv*U' ;

format long 

Fourier(1,:)

