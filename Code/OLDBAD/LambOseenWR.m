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

% W = 1;
% R = 1;
% Omega = 3.84/2;
% k = 6;%2*Omega/W;
% G = 5.9;
% d = 0.6;
% numPoints = 100;
% rstarguess = 0.5;
% %populate the params struct 
% params = struct();
% params.R = R;
% params.W = W;
% params.G = G;
% params.d = d;
params.k = k;
% params.Omega = Omega;   
% params.numPoints = numPoints;
r = linspace(0,R,numPoints)';
%100 pts from rstar - R

PsiRot = rotatingSol(1:end-1); %ignore the last term, rstar
rstar = rotatingSol(end);
rRotOld = [0, linspace(rstar,R)]';
PsiRotOld = [0;0;PsiRot(2:end)];
PsiRotguess = interp1(rRotOld,PsiRotOld,r);
epsilon = 1e-6;
PsiRotguess = PsiRotguess + epsilon * r.*(r-R);
%PsiRotguess = 0.5*r.^2;
% 
% PsiLam = lamSol(1:end-1); %ignore the last term, rstar
% rstarLam = lamSol(end);
% rRotOld = [0, linspace(rstarLam,R)]';
% PsiLamOld = [0;0;PsiLam(2:end)];
% PsiLamguess = interp1(rRotOld,PsiLamOld,r);
% 
% 
% [newlamSol,err,flag,struc] = fsolve(@(x)SLLambOseen(x,params),PsiLamguess,opt);
% figLO = PlotPsiAndW(lamSol,params);
% 

opt = optimoptions('fsolve','MaxFunctionEvaluations',1e6);
[newRotSol,err,flag,struc] = fsolve(@(x)SLRotating(x,params),PsiRotguess,opt);
figRot = PlotPsiAndW(newRotSol,params);

function fig = PlotPsiAndW(sol,params)
R = params.R;
numPoints=params.numPoints;

r = linspace(0,R,numPoints)';
Psi = sol;
% w = 1/r dpsidr
i = (2:numPoints-1)';
dr = r(2)-r(1);
dpsidr = (Psi(i+1) - Psi(i-1))/(2*dr);
w = dpsidr ./ r(i);
%w = (Psi(3:end) - Psi(1:end-2))./(2*r(2:end-1)*(r(2)-r(1)));
fig = figure;

subplot(2,1,1)
plot(r,Psi)
xlabel("r")
ylabel("$$\Psi$$")
axis([0,R,0,inf])
grid on

subplot(2,1,2)
plot(r(2:end-1),w)
xlabel("r")
ylabel("w")
axis([0,R,0,inf])
grid on


end

function dd = (y,x,linspaced)
if linspaced
   
end
    dPsidr = (y(3:end) - Psi(1:end-2))/(2*x(2:end) - x(1:end-1));
end
function res = SLRotating(Psi,params)
%Solve the SquireLong equation for rotating flow
%Squire Long equation in 1D is
%d2psidr2 - dpsidr/r = r^2 dHdPsi - C dCdPsi
%Inputs:
%x -> [Psi ; rstar]
%params -> Struct of parameters for the problem
%expects params to contain
%R, W, k, numPoints


%grab relevant parameters
R = params.R;
W = params.W;
k = params.k;
numPoints = params.numPoints;

%indexer for finite differences
i = (2:numPoints-1)';
r = linspace(0,R,numPoints)';
%linspace r -> dr is constant
dr = r(2)-r(1);
res = zeros(numPoints,1);
dPsidr = (Psi(i+1) - Psi(i-1))/(2*dr);
d2Psidr2 =  (Psi(i+1) - 2*Psi(i) + Psi(i-1))/(dr^2);

LHS =  d2Psidr2 - (dPsidr./(r(i)));
RHS =  (k^2/(2*W)) * r(i).^2 - k^2 .* Psi(i); 


RHS((Psi(i) < 0) | (Psi(i)>1/2)) = 0;
%wrstar = (Psi(2) - Psi(1))/dr;
% w = dPsidr./ r(i);
%residuals will be LHS-RHS for the middle points
%the psi(R) case, and the Psi(rstar) case
res(2:end-1) = LHS - RHS;
res(1) = 0;
res(end) = Psi(end) - 0.5*W*R^2;
%res(end) = wrstar;


end

function res = SLLambOseen(Psi,params)
%Solve the SquireLong equation for the Lamb Oseen Vortex
%Squire Long equation in 1D is
%d2psidr2 - dpsidr/r = r^2 dHdPsi - C dCdPsi


R = params.R;
W = params.W;
G = params.G;
d = params.d;
nPoints = length(Psi);
i = (2:nPoints-1)';
r = linspace(0,R,nPoints)';
dr = r(2)-r(1);
res = zeros(nPoints,1);
dPsidr = (Psi(i+1) - Psi(i-1))/(2*dr);
d2Psidr2 =  (Psi(i+1) - 2*Psi(i) + Psi(i-1))/(dr^2);
LHS =  d2Psidr2 - (dPsidr./r(i));
%There is a 0/0 limit as Psi -> 0
%Handle this by giving the taylor series 
%to the small values of Psi
RHSsmallPsi = ((G^2.*r.^2)./(2*W*d^4*pi^2) - (G^2.*Psi.*((2*((4*r.^2)/d^2 + 2))/(W*d^2) - (2*r.^2)/(W*d^4)))./(4*W*d^2*pi^2));
RHSsmallPsi(Psi > 0.0001| Psi <0) = 0;
RHSbigPsi = (G^2/(2*W*d^2*pi^2).* ((r.^2*W ./(2*Psi)) -1).*(exp(-2*Psi/(W*d^2)) - exp(-4*Psi/(W*d^2))));  
RHSbigPsi(Psi <=0.0001 | Psi<0) = 0;
RHS = RHSsmallPsi + RHSbigPsi;

%As per Rusak, Wang 
%if the RHS is negative - set it to 0
%RHS(RHS < 0) = 0;
%wrstar = (Psi(2) - Psi(1))/dr;
res(2:end-1) = LHS - RHS(i);
res(1) = Psi(1);
res(end) = Psi(end) - 0.5*W*R^2;
% res(end) = wrstar;
%%in case i want to see the full error
%dot(res,res)

end