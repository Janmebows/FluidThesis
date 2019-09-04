close all
R = 1;
W = 1;
G = 5.9;
d = 0.6;
rstarguess = 0.5;

%Make plots less repulsive
set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');

rstarFunc =@(rstarguess) rstarsolve(rstarguess,R,W,G,d);
%make it easier to call

[rstar,err,flag,struc] = fzero(rstarFunc,rstarguess);
if flag >0 && err ~= -0.5*W*R^2
    func = @(eta,Psi,rstar)psiODELambOseen(eta,Psi,rstar,R,W,G,d);
    [r,psi] = ode45(@(r,Psi)func(r,Psi,rstar),[rstar,R],[0,0]);

    w = psi(:,2)./r;
    plot(r,psi(:,1))
    ylabel("$$\Psi$$")
    axis([0,R,0,inf])
    figure
    plot(r,w)
    ylabel("$$w$$")
    axis([0,R,0,inf])
else
    fprintf('failed\n\n')
end

function out = rstarsolve(rstarguess,R,W,G,d)
if rstarguess <= 0 
    fprintf("Trying to leave bounds - bad fzero!\n")
    err = -10;
elseif rstarguess >= 1
    fprintf("Trying to leave bounds - bad fzero!\n")
    err = 10;
else
func = @(eta,Psi,rstar)psiODELambOseen(eta,Psi,rstar,R,W,G,d);
    [~,psi] = ode45(@(r,Psi)func(r,Psi,rstarguess),[rstarguess,R],[0,0]);
    err = psi(end,1) - 0.5*W*R^2;
end    
    out = err;

end
function out= psiODELambOseen(r,Psi,rstar,R,W,G,d)

  if Psi(1,:) < 0.0001
      rhs = (G^2.*r.^2)./(2*W*d^4*pi^2) - (G^2.*Psi(1,:).*((2*((4*r^2)/d^2 + 2))/(W*d^2) - (2*r.^2)/(W*d^4)))/(4*W*d^2*pi^2);
       %rhs = (G^2*r^2)/(2*W*d^4*pi^2);
  else
     rhs = G^2/(2*W*d^2*pi^2)* ((r.^2*W ./(2*Psi(1,:))) -1).*(exp(-2*Psi(1,:)/(W*d^2)) - exp(-4*Psi(1,:)/(W*d^2)));  
  end
dPsidr = Psi(2,:);
d2Psidr2 = dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end