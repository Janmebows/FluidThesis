close all
clear all
%Make plots less repulsive
set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');

step = 0.001;
R = 1;
r0= 0.99;
W = 1;
%guess = [rstar,rhat,k]
guess = [0.3 0.65 5];
options = optimoptions('fsolve','MaxFunctionEvaluations',1e5);
[sol,fval,exitflag]  = fsolve(@(x) wAndPsi(x,R,r0,W),guess,options);
workingVals=[];

while(r0 > 0)%&& isreal(sol))
    if(exitflag > 0 && isreal(sol))
        workingVals = [workingVals ;sol,r0];
        fval;
    end
    
    r0=r0-step;
    [sol,fval,exitflag]  = fsolve(@(x) wAndPsi(x,R,r0,W),guess,options);


% rstar = x(1);
% rhat = x(2);
% k =x(3);
end
%plot(r0,kr0)

plot(workingVals(:,4),workingVals(:,3).*workingVals(:,4))
xlabel('$$r_0$$')
ylabel('$$kr_0$$')
grid on

title("Rankine body solution space")
%mimic the axes from Escudier/Keller
axis([0,1,0,3.832])
saveas(gcf,'rankineSolutionSet.eps','epsc')
function out = wAndPsi(x,R,r0,W)
rstar = x(1);
rhat = x(2);
k =x(3);
Wh = W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 );
psi_s = -0.5*W*rstar.^2;
psi_h = 0.5*W*(r0.^2 - rhat.^2);
%Ad,Bc ->psi_inner_rstar = 0 ,
%      ->psi_inner_rhat = 0.5*W*r0^2
Ad =  rhat*bessely(1,k*rhat).*psi_s - rstar*bessely(1,k*rstar).*psi_h;
Bd = -rhat*besselj(1,k*rhat).*psi_s + rstar*besselj(1,k*rstar).*psi_h;
deter = rhat*rstar*(besselj(1,k*rstar)*bessely(1,k*rhat) - bessely(1,k*rstar)*besselj(1,k.*rhat));

%%%out(1) -> w(rstar)*det = 0
%%%out(2) -> w_outer(rhat) = w_inner(rhat)
%%%out(3) -> net momentum = 0
out(1) =  W + k*(Ad*besselj(0,k*rstar) + Bd*bessely(0,k*rstar))./deter;
out(2) =  W + k*(Ad*besselj(0,k*rhat)  + Bd*bessely(0,k*rhat) )./deter - Wh;
out(3) = -rhat^2 + 0.25*((rhat^4 - rstar^4)/(r0^2)) + (0.75* r0^2) + (0.5*r0^2*log(rhat^2/r0^2));
end

