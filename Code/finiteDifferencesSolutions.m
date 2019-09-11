%%Numerically solve the psi problem.
%psir2 - psir /r = 0

%psir = (psi(i+1) - psi(i-1))/(2*dr)
%psir2(i) = (psi(i+1) - 2*psi(i) + psi(i-1))/(dr^2)
clear all
close all

%%Normal Grid
R = 1;
W = 1;
rstar = 0.4;
nPoints = 3001;
r = linspace(rstar,R,nPoints)';
dr = r(2) - r(1);

%A is the coefficient matrix for psi
%Solve A*psi = b

i = 2:nPoints-1;

invdr = 1/dr;
A=...
sparse(i,i,-2*(invdr^2),nPoints,nPoints) ... 
+sparse(i,i-1,invdr^2 + 0.5*invdr./r(i),nPoints,nPoints) ...
+sparse(i,i+1,invdr^2 - 0.5*invdr./r(i),nPoints,nPoints) ...
+sparse([1 nPoints],[1 nPoints],1,nPoints,nPoints);
b = zeros(nPoints,1);
b(end) = 0.5*W*R^2;

psi = A\b;
plot(r,psi)
axis([0,1,-inf,inf])
w = zeros(size(psi));
w(i) = 0.5*invdr*(psi(i+1) - psi(i-1))./r(i);
hold on
plot(r,w,'r')

%%
%%homogeneous ODE for PSI
%d2Psi - dPsi/r + Psi(k^2 - 1/r^2) = 0
k=4.1647;
APsi=...
sparse(i,i,-2*invdr^2 +(k^2 - 1./(r(i).^2)),nPoints,nPoints) ... 
+sparse(i,i-1,(invdr^2 - 0.5*invdr)./r(i),nPoints,nPoints) ...
+sparse(i,i+1,(invdr^2 + 0.5*invdr)./r(i),nPoints,nPoints) ...
+sparse([1 nPoints],[1 nPoints],1,nPoints,nPoints);
%boundary conditions
rstarindex = 10;
%bPsi = zeros(nPoints,1);
bPsi = sparse(nPoints,1);
bPsi(1) = 0;
bPsi(end) = 0.5*W*R^2% + 0.5*W*R^3;


Psi = APsi\bPsi;
psi2 = Psi./r - 0.5*W*r;
w2 = zeros(size(psi2));
w2(i) = 0.5*invdr*(psi2(i+1) - psi2(i-1))./r(i);

plot(r,Psi);
%of course w = 1/r dpsi dr 
%%Staggered Grid


%%
