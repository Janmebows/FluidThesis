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


R = 1;
r0= 0.8;
W = 1;
%guess = [rstar,rhat,k]
guess = [0.1718 0.8686 3.8874];
% lb = [0,0,0];
% ub = [R,R,100];
%sol = lsqnonlin(@(x) wAndPsi(x,R,r0,W),guess,lb,ub);
[sol,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(@(x) wAndPsi(x,R,r0,W),guess);

    
    
rstar = sol(1)
rhat= sol(2)
k = sol(3)


r = linspace(0,R);
% C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
% Ad =  rstar.*bessely(1,k.*rstar).*(2*C -W) +k.*bessely(0,k.*rhat).*(0.5*W*rstar.^2);
% Bd = -rstar.*besselj(1,k.*rstar).*(2*C-W)-k.*besselj(0,k.*rhat).*(0.5*W*rstar.^2);
% deter = k.*rstar.*(besselj(0,k.*rhat).*bessely(1,k.*rstar) - bessely(0,k.*rhat).*besselj(1,k.*rstar));
% 
% wind = W*deter + k.*(Ad.*besselj(0,k.*r) + Bd.*bessely(0,k.*r));
% win = wind/deter;
% wout = 2*C*ones(size(r));

Wh = W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 );
psi_s = -0.5*W*rstar.^2;
psi_h = 0.5*W*(r0.^2 - rhat.^2);
Ad =  rhat*bessely(1,k*rhat).*psi_s - rstar*bessely(1,k*rstar).*psi_h;
Bd = -rhat*besselj(1,k*rhat).*psi_s + rstar*besselj(1,k*rstar).*psi_h;
deter = rhat*rstar*(besselj(1,k*rstar)*bessely(1,k*rhat) - bessely(1,k*rstar)*besselj(1,k.*rhat));
win =  W + k*(Ad*besselj(0,k*r) + Bd*bessely(0,k*r))./deter;
wout = Wh*ones(size(r));

wfull = [win(r<rhat), wout(r>=rhat)];
wfull(r<rstar) = 0;
% A = Ad./deter;
% B = Bd./deter;
% psiinner = 0.5*W*r.^2 + r.^(A.*besselj(1,k*r) + B.*bessely(1,k*r));
% D = 0.5*W*R^2 * (r0^2-rhat^2)/(R^2-rhat^2);
% psiouter = C.*r.^2 + D;
plot(r,wfull)

grid on
xlabel('$r$')
ylabel('$w$')

psiinner = 0.5*W*r.^2 + r.*(Ad.*besselj(1,k*r) + Bd.*bessely(1,k*r))./deter;
psiouter = 0.5*Wh*r.^2 + 0.5*(W - Wh)*R^2;
psifull = [psiinner(r<rhat), psiouter(r>=rhat)];
psifull(r<rstar) = 0;
figure
plot(r,psifull)

grid on
xlabel('$r$')
ylabel('$\psi$')

kr0 = k*r0


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

