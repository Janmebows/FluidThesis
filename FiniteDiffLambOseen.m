%%solve lamb oseen with finite diff + fsolve
%%Numerically solve the psi problem.
%psir2 - psir /r = f(psi,r)
clear all
close all
%Make plots less repulsive
set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');

%%Normal Grid
R = 1;
W = 1;
rstar = 0.3;
nPoints = 101;
r = linspace(rstar,R,nPoints)';
dr = r(2) - r(1);
%eta = linspace(0,1,nPoints)';
%deta = eta(2)-eta(1);

%dpsidr = (psi(i+1) - psi(i-1))/(2*dr)
%d2psidr2(i) = (psi(i+1) - 2*psi(i) + psi(i-1))/(dr^2)




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

%options=optimset('FinDiffType','central','MaxIter',10000);
%Sol = fsolve(, options);



psi = A\b;
plot(r,psi)
axis([0,1,-inf,inf])
w = zeros(size(psi));
w(i) = 0.5*invdr*(psi(i+1) - psi(i-1))./r(i);


function out = lambOseen(Psi)
f =@(r,psi) A * ((r.^2*B./(2*Psi(1,:)) -1).*(exp(-2*Psi(1,:)*C) - exp(-4*Psi(1,:)*C)));
g =@(r,psi) 2*dr^2 * f(r,psi) + 4*psi;

end


