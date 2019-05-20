%close all
%clear all
W = 1;
R = 1;
rstar = 0.9;
kguess = 1;
krstar = [];
[k,fval,exitflag,output,jacobian] =fsolve(@(k)wrstar(k,rstar,W,R),kguess,optimoptions('fsolve','OptimalityTolerance', 1e-5)); 

while(exitflag > 0 && isreal(k))
krstar = [krstar;k,rstar];
rstar = rstar - 0.01;
[k,fval,exitflag,output,jacobian] =fsolve(@(k)wrstar(k,rstar,W,R),kguess,optimoptions('fsolve','OptimalityTolerance', 1e-5)); 
end
figure
plot(krstar(:,1),krstar(:,2))
xlabel('k')
ylabel('r_*')
axis([0,10,0,1])
title("Numerically obtained solution set")
saveas(gcf,'simpleProblemNumericalSolutionSpace.eps','epsc')
function out = wrstar(k,r,W,R)
deter = @(r,k,R)(besselj(1,k.*r).*bessely(1,k.*R) - besselj(1,k.*R).*bessely(1,k.*r));
Atimesdet =@(r,k,R) -0.5*W*r .*bessely(1,k*R) ;
Btimesdet =@(r,k,R) 0.5*W*r.*besselj(1,k*R);
out = W*deter(r,k,R) + k.*(Atimesdet(r,k,R).*besselj(0,k.*r) + Btimesdet(r,k,R).*bessely(0,k.*r));
end


