function w = Plotw(params)
%Make plots less repulsive
set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');

rstar = params(1);
rhat = params(2);
k = params(3);
r0=params(4);
R=params(5);
W= params(6);
r = linspace(0,R);

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
plot(r,wfull)
hold on
plot(rhat,1,'x')
plot(rstar,0,'x')
end
