%%Numerically solve the psi problem.
%psir2 - psir /r = 0
close all
clear all
%psir = (psi(i+1) - psi(i-1))/(2*dr)
%psir2(i) = (psi(i+1) - 2*psi(i) + psi(i-1))/(dr^2)
figure

%%Normal Grid
R = 1;
W = 1;
%rstar = 0.01;
nPoints = 501;
r = linspace(0,R,nPoints)';
dr = r(2) - r(1);
%special k
k = 3.831706186854435;
Omega = k*W/2;
%A is the coefficient matrix for psi
%Solve A*psi = b

i = 2:nPoints-1;

invdr = 1/dr;


A=...
sparse(i,i,-2*(invdr^2) + k^2,nPoints,nPoints) ... 
+sparse(i,i-1,invdr^2 + 0.5*invdr./r(i),nPoints,nPoints) ...
+sparse(i,i+1,invdr^2 - 0.5*invdr./r(i),nPoints,nPoints) ...
+sparse([1 nPoints],[1 nPoints],1,nPoints,nPoints);
b = zeros(nPoints,1);
b(i) = 2*r(i).^2 * Omega^2/W;
b(end) = 0.5*W*R^2;

psi = A\b;
plot(r,psi)

w = zeros(size(psi));
w(i) = 0.5*invdr*(psi(i+1) - psi(i-1))./r(i);
%have to do 1 and end "manually"
w(1) = 2*invdr^2 *(psi(2) - psi(1));
%w(end) =  
plot(r,w)
grid on