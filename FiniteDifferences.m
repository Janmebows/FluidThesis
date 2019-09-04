close all
%Make plots less repulsive
set(groot, 'DefaultLineLineWidth', 1, ...
    'DefaultAxesLineWidth', 1, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultColorbarTickLabelInterpreter', 'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');

W = 1;
R = 1;
Omega = 2.5;
k = 2*Omega/W;
G = 1;
d = 1;
rstarguess = 0.5;
params = struct();
params.R = R;
params.W = W;
params.G = G;
params.d = d;
params.k = k;
params.Omega = Omega;
func = @rotatingEquation;
guessFunc = @(rstar)findrstar(func, rstar,params);
[rstar,err,flag,struc] = fzero(guessFunc,[0,1]);
[r,Psi,~] = rotatingEquation(rstar,params);
i = 2:length(Psi)-1;
w = (1./r(i)).*(Psi(i+1) - Psi(i-1))/(2*(r(2)-r(1)));

plot(r,Psi)
axis([0,1,0,inf])
figure
plot(r(i),w)
axis([0,1,0,inf])


function out = findrstar(func, rstar,params)
if rstar <= 0 
    fprintf("Trying to leave bounds - too low - bad fzero!\n")
    err = -10;
elseif rstar >= 1
    fprintf("Trying to leave bounds - too high - bad fzero!\n")
    err = 10;
else
    [~,~,err] = func(rstar,params);
    err = dot(err,err);
end
    out = err;
    

end
function [r,Psi,err] = rotatingEquation(rstar,params)
R = params.R;
W = params.W;
%G = params.G;
%d = params.d;
Omega = params.Omega;
k = params.k;
nPoints = 100;
i = (2:nPoints-1)';
%eta = linspace(0,1,nPoints)';
%r = eta*(R-rstar)+rstar;
r = linspace(rstar,R,nPoints)';
dr = r(2)-r(1);%eta(2)-eta(1);

%solves the system with eta instead of r
%ROTATING FLOW
%REPLACE ALL r TERMS WITH eta TERMS
dPsidr = @(Psi) (Psi(i+1) - Psi(i-1))/(2*dr);
d2Psidr2 = @(Psi) (Psi(i+1) - 2*Psi(i) + Psi(i-1))/(dr^2);
LHS = @(Psi) d2Psidr2(Psi) - (dPsidr(Psi)./r(i));
%RHS = @(Psi) 2*r(i).^2*Omega^2/W - 4*Psi(i)*Omega^2/W^2;
RHS = @(Psi) (k^2/(2*W)) * r(i).^2 - k^2 .* Psi(i); 
%RHS = @(Psi) (k^2/(2*W)).* r(i).^2 - (k^2 .* Psi(i));
%Equations:
%Psi(1) =0
%Psi(end) = 0.5*W*R^2
%Psi(2:end-1) = ...
eqn = @(Psi) [Psi(1); LHS(Psi)-RHS(Psi); Psi(end)-0.5*W*R^2];
initGuess = 0.5*W*r.^2;

options = optimset('Display','off');
[Psi,err,flag,struc] = fsolve(eqn,initGuess,options);

end