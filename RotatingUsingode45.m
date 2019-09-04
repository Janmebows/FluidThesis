close all
R = 1;
W = 1;
k = 3.84;
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



func = @(r,Psi,rstar)psiODE(r,Psi,rstar,R,W,k);

rstarFunc =@(rstarguess) rstarsolve(rstarguess,R,W,k);
%make it easier to call

[rstar,err,flag,struc] = fzero(rstarFunc,rstarguess);

[r,Psi] = ode45(@(r,Psi)func(r,Psi,rstar),[rstar,R],[0,0]);


w = Psi(:,2)./r;
plot(r,Psi(:,1))
ylabel("$$\Psi$$")
axis([0,R,0,inf])
figure
plot(r,w)
ylabel("$$w$$")
axis([0,R,0,inf])




function out = rstarsolve(rstarguess,R,W,k)
if rstarguess <= 0 
    fprintf("Trying to leave bounds - bad fzero!\n")
    err = -10;
elseif rstarguess >= 1
    fprintf("Trying to leave bounds - bad fzero!\n")
    err = 10;
else
    func = @(r,Psi)psiODE(r,Psi,rstarguess,R,W,k);
    [~,psi] = ode45(@(r,Psi)func(r,Psi),[rstarguess,R],[0,0]);

    err = psi(end,1) - 0.5*W*R^2;
end    
    out = err;

end
function out= psiODE(r,Psi,rstar,R,W,k)
rhs = (k^2/(2*W)).* r.^2 - k^2 .* Psi(1,:); 
dPsidr = Psi(2,:);
d2Psidr2 =dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end