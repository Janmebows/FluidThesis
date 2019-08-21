syms u v w r G d C H Psi(r) R W 
u = 0;
v = G/(2*pi*r) * (1-exp(-r^2/d^2))
psi = 0.5*W*r^2;
w = W
[psiEq,rhs] = SymSquireLong1D(u,v,w,psi)

dPsi = diff(Psi(r))
cond = [Psi(0) == 0, Psi(R) == 0, dPsi(0) ==0];
dsolve(psiEq,Psi,cond)
%dsolve(psiEq)


%%
close all
R = 1;
W = 1;
G = 3;
d = 0.05;
rstarguess = 0.1;

%Make plots less repulsive
set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');


%singular de must include the singular point
%in this case dpsi/r(0) = 0
initguess = @(eta) [0.5*W*eta^2,W*eta]%[1.5*eta.^2+1,1.5*eta.^2];
%may need to find a way to solve right to left rather than left to right.
numPts = 100;
solinit=bvpinit(linspace(0,1,numPts),initguess,rstarguess);
options = []%bvpset('RelTol',1e-5,'AbsTol',1e-8,'Nmax',1e6,'Vectorized','on');%,'SingularTerm',[0,0;0,1]);
func = @(eta,Psi,rstar)psiODE(eta,Psi,rstar,R,W,G,d);
bound = @(psirstar,psiR,rstar)boundaries(psirstar,psiR,rstar,R,W);
sol=bvp5c(func,bound,solinit,options);
rstar = sol.parameters;
eta = sol.x;
r = etaTor(sol.x,rstar,R);
psi = sol.y(1,:);
w = sol.y(2,:)./r;
plot(r,psi)
ylabel("$$\Psi$$")
figure
plot(r,w)
ylabel("$$w$$")




function res = boundaries(psirstar,psiR,rstar,R,W)
%the boundary conditions
%psi(r*) = 0
%psi(R) = 0.5*W*R^2
%w(r*) = dPsidr/r(r*) =0
res = [psirstar(1);
        psiR(1) - 0.5*W*R^2;
        psirstar(2)];
end

function r = etaTor(eta,rstar,R)
%eta is 0,...,1
%convert to r = rstar, ..., R
r = eta.*(R-rstar) + rstar;

end

function eta = rToeta(r,rstar,R)
eta = (r-rstar)./(R-rstar);
end

function out= psiODE(eta,Psi,rstar,R,W,G,d)
%dpsi2dr2 - dpsidr / r = f(r,psi)


%at the starting point, Psi = dPsi = 0 and so rhs = NaN. Have to deal with
%the division by Psi.
%convert eta to r
%eta = (r-rstar)/(R-rstar)

r = etaTor(eta,rstar,R);

%  if Psi(1,:) < 0.0001
%     PsiSeries = 1 - (-1+Psi(1,:)) + (-1+Psi(1,:)).^2;
%     rhs = G^2/(2*W*d^2*pi^2) * ((r.^2*W .*PsiSeries) -1).*(exp(-2*Psi(1,:)/(W*d^2)) - exp(-4*Psi(1,:)/(W*d^2))); 
%  else
%     rhs = G^2/(2*W*d^2*pi^2) * ((r.^2*W./(2*Psi(1,:)) -1).*(exp(-2*Psi(1,:)/(W*d^2)) - exp(-4*Psi(1,:)/(W*d^2))));
% 
%  end
%homo case
rhs = (2*G^2/W)* r^2 - (4*G^2/W^2) * Psi(1,:); 
%rhs = Psi(1,:) .*(-4*G^2/W^2 + 1./r.^2);

dPsidr = Psi(2,:);

dPsi2dr2 = dPsidr./r + rhs;

%out = [dPsidr ; dPsi2dr2 ; w];
%out = [dPsidr ; dPsi2dr2;dPsidr/r];
out = [dPsidr ; dPsi2dr2];
end