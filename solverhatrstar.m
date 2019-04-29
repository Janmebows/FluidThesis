close all
clear all
R = 1;
r0= 0.5;
rstarguess=0.1;
W = 1;
k = 20;
choosek = false;
%Guess takes an initial guess for [rhat,rstar,k]
guess = [0.5,rstarguess,k];
%prob.objective = @(x) wAndPsi([x,k],R,r0,W);
%prob.x0 = guess(1:2);
%prob.lb = [0,0];
%prob.ub = [R,R];
%prob.solver = 'fsolve';
%prob.options = [];
%fsolve(prob)
if(choosek)
    guess(3)=[];
    [sol,val,exitflag,output,jacobian] = fsolve(@(x) wAndPsi([x,k],R,r0,W),guess);
else
    [sol,val,exitflag,output,jacobian] = fsolve(@(x) wAndPsi(x,R,r0,W),guess);
end
if(complex(sol(1))~=0 || complex(sol(2)) ~=0)
%    matException= MException('VerifyOutput:OutOfBounds',"Result cannot be complex\n try a different input");
%    throw(matException)
end
exitflag
rhat = sol(1);
rstar = sol(2);

r = linspace(0,R);
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad =  rstar.*bessely(1,k.*rstar).*(2*C -W) -k.*bessely(0,k.*rhat).*(0.5*W*rstar.^2);
Bd = -rstar.*besselj(1,k.*rstar).*(2*C-W)-k.*besselj(0,k.*rhat).*(0.5*W*rstar.^2);
deter = k.*rstar.*(besselj(0,k.*rhat).*bessely(1,k.*rstar) - bessely(0,k.*rhat).*besselj(1,k.*rstar));

wind = W*deter + k.*(Ad.*besselj(0,k.*r) + Bd.*bessely(0,k.*r));
win = wind/deter;
wout = 2*C*ones(size(r));
wfull = [win(r<rhat), wout(r>=rhat)];
plot(r,wfull)


function out = wAndPsi(x,R,r0,W)
rhat = x(1);
rstar = x(2);
k = x(3);
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad =  rstar.*bessely(1,k.*rstar).*(2*C-W)-k.*bessely(0,k.*rhat).*(0.5*W*rstar.^2);
Bd = -rstar.*besselj(1,k.*rstar).*(2*C-W)-k.*besselj(0,k.*rhat).*(0.5*W*rstar.^2);
deter = k.*rstar.*(besselj(0,k.*rhat).*bessely(1,k.*rstar) - bessely(0,k.*rhat).*besselj(1,k.*rstar));

%%%out(1) -> w(rstar)*det = 0
%%%out(2) -> psi(rhat)*det = 0.5*W*r0^2*det
out(1) =  W*deter + k.*(Ad.*besselj(0,k.*rstar) + Bd.*bessely(0,k.*rstar));
out(2) = 0.5*W*deter.*(rhat.^2-r0^2) + rhat.*(Ad.*besselj(1,k.*rhat) + Bd.*bessely(1,k.*rhat))- (0.5*deter*W*r0^2);

end

