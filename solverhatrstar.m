close all
clear all
R = 1;
r0= 0.9;
W = 1;
k = 10;
choosek = false;
guess = [r0,0.1,10];
if(choosek)
    guess(3)=[];
    [sol,val] = fsolve(@(r) wAndPsi([r,k],R,r0,W),guess);
else
    [sol,val] = fsolve(@(r) wAndPsi(r,R,r0,W),guess);

end

rhat = sol(1);
rstar = sol(2);

r = linspace(0,R);
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad = (rstar.*bessely(1,k.*rstar)).*(2*C -W) -k.*bessely(0,k.*rhat).*(0.5*W*rstar.^2);
Bd = -rstar.*besselj(1,k.*rstar).*(2*C-W)+k.*besselj(0,k.*rhat).*(0.5*W*rstar.^2);
deter = k.*rstar.*(besselj(0,k.*rhat).*bessely(1,k.*rstar) - bessely(0,k.*rhat).*besselj(1,k.*rstar));

wind = W*deter + k.*(Ad.*besselj(0,k.*r) + Bd.*bessely(0,k.*r));
win = wind/deter;
wout = 2*C*ones(size(r));
wfull = [win(r<rhat), wout(r>=rhat)];
plot(r,wfull)
function out = wAndPsi(r,R,r0,W)
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - r(1).^2 ) ;
Ad = (r(2).*bessely(1,r(3).*r(2))).*(2*C -W) -r(3).*bessely(0,r(3).*r(1)).*(0.5*W*r(2).^2);
Bd = -r(2).*besselj(1,r(3).*r(2)).*(2*C-W)+r(3).*besselj(0,r(3).*r(1)).*(0.5*W*r(2).^2);
deter = r(3).*r(2).*(besselj(0,r(3).*r(1)).*bessely(1,r(3).*r(2)) - bessely(0,r(3).*r(1)).*besselj(1,r(3).*r(2)));

out(1) =  W*deter + r(3).*(Ad.*besselj(0,r(3).*r(2)) + Bd.*bessely(0,r(3).*r(2)));
out(2) = 0.5*W*deter.*(r(1).^2-r0^2) + r(1).*(Ad.*besselj(1,r(3).*r(1)) + Bd.*bessely(1,r(3).*r(1)))- (0.5*deter*W*r0^2);

end

