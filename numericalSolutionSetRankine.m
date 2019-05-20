close all
clear all
R = 1;
r0= 0.3;
rstar=0.01;
W = 1;
guess = [0.6,10];
options = optimoptions('fsolve','MaxFunctionEvaluations',1e3);
[sol,fval,exitflag]  = fsolve(@(x) wAndPsi(x,rstar,R,r0,W),guess,options);
workingVals=[];

while(exitflag > 0 && isreal(sol))
workingVals = [workingVals ;sol(2),rstar,sol(1)];
[sol,fval,exitflag]  = fsolve(@(x) wAndPsi(x,rstar,R,r0,W),guess,options);
rstar=rstar+0.01;

%rhat= sol(1);
%k= sol(2);
end
plot(workingVals(:,1),workingVals(:,2))



function out = wAndPsi(x,rstar,R,r0,W)
%rstar = x(1);
rhat = x(1);
k = x(2);
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad =  rstar*bessely(1,k*rstar)*(2*C-W)+k*bessely(0,k*rhat)*(0.5*W*rstar^2);
Bd = -rstar*besselj(1,k*rstar).*(2*C-W)-k*besselj(0,k*rhat).*(0.5*W*rstar^2);
deter = k*rstar*(besselj(0,k*rhat)*bessely(1,k*rstar) - bessely(0,k*rhat)*besselj(1,k.*rstar));

%%%out(1) -> w(rstar)*det = 0
%%%out(2) -> psi(rhat)*det = 0.5*W*r0^2*det
out(1) =  W*deter + k*(Ad*besselj(0,k*rstar) + Bd*bessely(0,k*rstar));
out(2) = 0.5*W*deter*(rhat.^2-r0^2) + rhat*(Ad.*besselj(1,k*rhat) + Bd*bessely(1,k*rhat))- (0.5*deter*W*r0^2);

end

