close all
clear all
R = 1;
r0= 0.5;
rstarguess=0.2;
W = 1;
guess = [rstarguess,0.6,5*sqrt(2)];
lb = [0,0,0];
ub = [R,R,100];
%sol = lsqnonlin(@(x) wAndPsi(x,R,r0,W),guess,lb,ub);
[sol,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(@(x) wAndPsi(x,R,r0,W),guess);

    
rstar = sol(1);
rhat= sol(2);
k = sol(3);
r = linspace(0,R);
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad =  rstar.*bessely(1,k.*rstar).*(2*C -W) +k.*bessely(0,k.*rhat).*(0.5*W*rstar.^2);
Bd = -rstar.*besselj(1,k.*rstar).*(2*C-W)-k.*besselj(0,k.*rhat).*(0.5*W*rstar.^2);
deter = k.*rstar.*(besselj(0,k.*rhat).*bessely(1,k.*rstar) - bessely(0,k.*rhat).*besselj(1,k.*rstar));

wind = W*deter + k.*(Ad.*besselj(0,k.*r) + Bd.*bessely(0,k.*r));
win = wind/deter;
wout = 2*C*ones(size(r));
wfull = [win(r<rhat), wout(r>=rhat)];
wfull(r<rstar) = 0;
A = Ad./deter;
B = Bd./deter;
psiinner = 0.5*W*r.^2 + r.^(A.*besselj(1,k*r) + B.*bessely(1,k*r));
D = 0.5*W*R^2 * (r0^2-rhat^2)/(R^2-rhat^2);
psiouter = C.*r.^2 + D;
psifull = [psiinner(r<rhat), psiouter(r>=rhat)];
plot(r,wfull)
%plot(r,psifull)

function out = wAndPsi(x,R,r0,W)
rstar = x(1);
rhat = x(2);
k =x(3);
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad =  rstar*bessely(1,k*rstar)*(2*C-W)+k*bessely(0,k*rhat)*(0.5*W*rstar^2);
Bd = -rstar*besselj(1,k*rstar).*(2*C-W)-k*besselj(0,k*rhat).*(0.5*W*rstar^2);
deter = k*rstar*(besselj(0,k*rhat)*bessely(1,k*rstar) - bessely(0,k*rhat)*besselj(1,k.*rstar));

%%%out(1) -> w(rstar)*det = 0
%%%out(2) -> psi(rhat)*det = 0.5*W*r0^2*det
%%%out(3) -> net momentum = 0
out(1) =  W*deter + k*(Ad*besselj(0,k*rstar) + Bd*bessely(0,k*rstar));
out(2) = 0.5*W*deter*(rhat.^2-r0^2) + rhat*(Ad.*besselj(1,k*rhat) + Bd*bessely(1,k*rhat))- (0.5*deter*W*r0^2);
out(3) = -rhat^2 + 0.25*((rhat^4 - rstar^4)/(r0^2)) + (0.75* r0^2) + (0.5*r0^2*log(rhat^2/r0^2));
end

