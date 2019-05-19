close all 
clear all
global W k R r0
W = 1;
k = 5;
R = 1;
r0 = 0.50;

det =@(rhat,rstar) rhat*(bessely(1,k*rhat)*besselj(1,k*rstar) - besselj(1,k*rhat)*bessely(1,k*rstar));
A =@(rhat,rstar) 1/det(rhat,rstar) * (rhat*bessely(1,k*rhat)*(-0.5*W*rstar) -bessely(1,k*rstar)*(0.5*W*(r0^2-rhat^2)));
B =@(rhat,rstar) 1/det(rhat,rstar) * (-rhat*besselj(1,k*rhat)*(-0.5*W*rstar) +besselj(1,k*rstar)*(0.5*W*(r0^2-rhat^2)));
C =@(rhat,rstar) 0.5*(W + k*(A(rhat,rstar)*besselj(0,k*rhat) + B(rhat,rstar)*bessely(0,k*rhat)));
D =@(rhat,rstar) 0.5*W*R^2 -C(rhat,rstar);
rhateqn =@(rhat,rstar) (A(rhat,rstar)*besselj(0,k*rhat) + B(rhat,rstar)*bessely(0,k*rhat))*(rhat^2 -1) == W/k * (r0^2 - R^2 -1);
rstareqn =@(rhat,rstar) W + k*(A(rhat,rstar)*besselj(0,k*rstar) + B(rhat,rstar)*bessely(0,k*rstar));

x0 = [0.5,0.3];
lb = [0,0];
ub = [R,R];
[x,resnorm,residual,flag] = lsqnonlin(@equations,x0,lb,ub);

rhat = x(1);
rstar = x(2);
r= linspace(0,R);
PsiIn = 0.5*W*rstar.^2 + r.*(A(rhat,rstar) .* besselj(1,k*r) + B(rhat,rstar).*bessely(1,k*r));
PsiOut = C(rhat,rstar)* r.^2 + D(rhat,rstar);
Psifull = [PsiIn(r<rhat),PsiOut(r>=rhat)];
%plot(r,Psifull)

Winfunc =@(r) W + k*(A(rhat,rstar)*besselj(0,k*r) + B(rhat,rstar)*bessely(0,k*r));
Win = Winfunc(r);
Wout = ones(1,length(r))*2*C(rhat,rstar);

Wfull = [Win(r<rhat), Wout(r>=rhat)];
%Wfull(r<rstar) = 0;
Wfull(r<rstar) = Winfunc(rstar);
figure
plot(r,Wfull)
function output = equations(x)
global W k R r0
rhat = x(1);
rstar = x(2);
det = rhat*(bessely(1,k*rhat)*besselj(1,k*rstar) - besselj(1,k*rhat)*bessely(1,k*rstar));
A = 1/det * (rhat*bessely(1,k*rhat)*(-0.5*W*rstar) -bessely(1,k*rstar)*(0.5*W*(r0^2-rhat^2)));
B = 1/det * (-rhat*besselj(1,k*rhat)*(-0.5*W*rstar) +besselj(1,k*rstar)*(0.5*W*(r0^2-rhat^2)));
C = 0.5*(W + k*(A*besselj(0,k*rhat) + B*bessely(0,k*rhat)));
D = 0.5*W*R^2 -C;
rhateqnval =(A*besselj(0,k*rhat) + B*bessely(0,k*rhat))*(rhat^2 -1) -W/k * (r0^2 - R^2 -1);
rstareqnval = W + k*(A*besselj(0,k*rstar) + B*bessely(0,k*rstar));
output = [rhateqnval,rstareqnval];

end