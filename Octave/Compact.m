p   = 5 ;   % band on lhs
N   = 8 ;   % Order of accuracy
pts = 9 ;

%  4th compact symmetric:  al = 1/4,              ; a = 3/4
%  6th compact Lele     :  al = 1/3,              ; a = 28/36,     b = 1/36
%  8th compact symmetric:  al = 16/36,  be = 1/36 ; a = 400/3/180, b = 1/180

s1 = (p-1)/(p+1) ;
s2 = (p-3)/(p+3) ;
s3 = (p-5)/(p+5) ;
s4 = (p-7)/(p+7) ;
t1 = (N - (p-1)) / (N - (p-3)) ;
t2 = (N - (p+1)) / (N - (p-5)) ;
t3 = (N - (p+3)) / (N - (p-7)) ;
t4 = (N - (p+5)) / (N - (p-9)) ;

al = s1*t1                   ;
be = s1*t1*s2*t2             ;
ga = s1*t1*s2*t2*s3*t3       ;
de = s1*t1*s2*t2*s3*t3*s4*t4 ;

r = (p-1)/2 ;
phi = 1 ;
for i = 1:r 
 phi = phi * (i/(4*i+2))
end
ga1 = 1 ;
for i = 1:(1-r)
  ga1 = ga1 * (-i/(p+i))
end
ga2 = 1 ;
for i = 1:(2-r)
  ga2 = ga2 * (-i/(p+i))
end
ga3 = 1 ;
for i = 1:(3-r)
  ga3 = ga3 * (-i/(p+i))
end

a = 4 / (p+1)^2 * (p*N - (p-1)^2)*(N+2) / (N-(p-3))^2 ;
b = phi*ga1*t1*t2 ;
c = phi*ga2*t1*t2*t3 ;
d = phi*ga3*t1*t2*t3*t4 ;

pts = 9 ;

qcol9 = [0,-a,-b,-c,-d,+d,+c,+b,+a] ;
qrow9 = [0,+a,+b,+c,+d,-d,-c,-b,-a] ;
pcol9 = [1,al,be,ga,de,de,ga,be,al] ;
qmat9 = toeplitz(qcol9,qrow9) ;
pmat9 = toeplitz(pcol9) ;
dmat9 = inv(pmat9)*qmat9/(2*pi/pts) ;
%  4th order compact
 qcolI = [0,-1,0,0,0,0,0,0,+1]/2;
 qrowI = [0,+1,0,0,0,0,0,0,-1]/2;
 pcolI = [4,+1,0,0,0,0,0,0,+1]/6;
%  6th order compact
 qcolI = [0,-28,-1,0,0,0,0,+1,+28]/36;
 qrowI = [0,+28,+1,0,0,0,0,-1,-28]/36;
 pcolI = [3,+1,0,0,0,0,0,0,+1]/3;
%  8th order compact
 qcolI = [0,-400/3,-1,0,0,0,0,+1,+400/3]/180;
 qrowI = [0,+400/3,+1,0,0,0,0,-1,-400/3]/180;
 pcolI = [36,+16,1,0,0,0,0,1,+16]/36;

qmatI = toeplitz(qcolI,qrowI) ;
pmatI = toeplitz(pcolI) ;
dmatI= inv(pmatI)*qmatI/(2*pi/pts) ;

%qcol2 = [0,-8,+1,0,0,0,0,-1,+8]/12;
%qrow2 = [0,+8,-1,0,0,0,0,+1,-8]/12;
%qcol2 = [0,-45,+9,-1,0,0,+1,-9,+45]/60;
%qrow2 = [0,+45,-9,+1,0,0,-1,+9,-45]/60;
qcolE = [0,-224,+56,-32/3,+1,-1,+32/3,-56,+224]/280;
qrowE = [0,+224,-56,+32/3,-1,+1,-32/3,+56,-224]/280;
dmatE = toeplitz(qcolE,qrowE)/(2*pi/pts) ;

%  Fourier discretization
qrowF = [0,+1.46190220008155e+00,-7.77861913430207e-01,+5.77350269189626e-01,-5.07713305942873e-01,+5.07713305942873e-01,-5.77350269189626e-01,+7.77861913430207e-01,-1.46190220008155e+00];
qcolF = [0,-1.46190220008155e+00,+7.77861913430207e-01,-5.77350269189626e-01,+5.07713305942873e-01,-5.07713305942873e-01,+5.77350269189626e-01,-7.77861913430207e-01,+1.46190220008155e+00];
dmatF = toeplitz(qcolF,qrowF) %/(2*pi/pts) ;

kk = 9 ;
vecs = zeros(pts,kk) ;
vecc = zeros(pts,kk) ;
for j = 1:kk
   for i = 1:pts
     x = 2*pi*(i-1)/(pts+0) ;
     vecs(i,j) =   sin(j*x) ;
     vecc(i,j) =   cos(j*x) ;
   end
end
tmps = dmatE*vecs ;
tmpc = dmatE*vecc ;
tmps = dmatI*vecs ;
tmpc = dmatI*vecc ;
tmps = dmatF*vecs ;
tmpc = dmatF*vecc ;

for j = 1:kk
  wrks(:,j) = tmps(:,j) - j*vecc(:,j) ;
  wrkc(:,j) = tmpc(:,j) + j*vecs(:,j) ;
end
