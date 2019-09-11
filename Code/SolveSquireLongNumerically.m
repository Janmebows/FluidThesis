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
Omega = 3.84/2;
k = 4.1647;%2*Omega/W;
G = 5.9;
d = 0.6;
numPoints = 100;
rstarguess = 0.5;
%populate the params struct 
params = struct();
params.R = R;
params.W = W;
params.G = G;
params.d = d;
params.k = k;
params.Omega = Omega;   
params.numPoints = numPoints;
r = linspace(0,R,numPoints)';
Psiguess = 0.5*r.^2;


[lamSol,err,flag,struc] = fsolve(@(x)SLLambOseen(x,params),[Psiguess; rstarguess]);
figLO = PlotPsiAndW(lamSol,params);
[rotatingSol,err,flag,struc] = fsolve(@(x)SLRotating(x,params),[Psiguess; rstarguess]);
figRot = PlotPsiAndW(rotatingSol,params);

saveas(figLO,'../Plots/LOSolFD.eps','epsc')
saveas(figRot,'../Plots/RotSolFD.eps','epsc')
function fig = PlotPsiAndW(sol,params)
R = params.R;
numPoints=params.numPoints;

rstar = sol(end);
r = linspace(rstar,R,numPoints)';
Psi = sol(1:end-1);
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

function res = SLRotating(x,params)
%Solve the SquireLong equation for rotating flow
%Squire Long equation in 1D is
%d2psidr2 - dpsidr/r = r^2 dHdPsi - C dCdPsi
%Inputs:
%x -> [Psi ; rstar]
%params -> Struct of parameters for the problem
%expects params to contain
%R, W, k, numPoints
%separate x into psi and rstar
Psi = x(1:end-1);
rstar = x(end);

%grab relevant parameters
R = params.R;
W = params.W;
k = params.k;
numPoints = params.numPoints;

%indexer for finite differences
i = (2:numPoints-1)';
r = linspace(rstar,R,numPoints)';
%linspace r -> dr is constant
dr = r(2)-r(1);
res = zeros(numPoints+1,1);
dPsidr = (Psi(i+1) - Psi(i-1))/(2*dr);
d2Psidr2 =  (Psi(i+1) - 2*Psi(i) + Psi(i-1))/(dr^2);
LHS =  d2Psidr2 - (dPsidr./r(i));
RHS =  (k^2/(2*W)) * r(i).^2 - k^2 .* Psi(i); 

RHS(RHS < 0) = 0;
wrstar = (Psi(2) - Psi(1))/dr;
res(2:end-2) = LHS - RHS;
res(1) = Psi(1);
res(end-1) = Psi(end) - 0.5*W*R^2;
res(end) = wrstar;


end
function res = SLLambOseen(x,params)
%Solve the SquireLong equation for the Lamb Oseen Vortex
%Squire Long equation in 1D is
%d2psidr2 - dpsidr/r = r^2 dHdPsi - C dCdPsi
Psi = x(1:end-1);
rstar = x(end);

R = params.R;
W = params.W;
G = params.G;
d = params.d;
nPoints = length(Psi);
i = (2:nPoints-1)';
r = linspace(rstar,R,nPoints)';
dr = r(2)-r(1);
res = zeros(nPoints+1,1);
dPsidr = (Psi(i+1) - Psi(i-1))/(2*dr);
d2Psidr2 =  (Psi(i+1) - 2*Psi(i) + Psi(i-1))/(dr^2);
LHS =  d2Psidr2 - (dPsidr./r(i));
%There is a 0/0 limit as Psi -> 0
%Handle this by giving the taylor series 
%to the small values of Psi
RHSsmallPsi = (Psi(i)<0.0001) .*((G^2.*r(i).^2)./(2*W*d^4*pi^2) - (G^2.*Psi(i).*((2*((4*r(i).^2)/d^2 + 2))/(W*d^2) - (2*r(i).^2)/(W*d^4)))./(4*W*d^2*pi^2));
RHSbigPsi = (Psi(i) >=0.0001) .* (G^2/(2*W*d^2*pi^2).* ((r(i).^2*W ./(2*Psi(i))) -1).*(exp(-2*Psi(i)/(W*d^2)) - exp(-4*Psi(i)/(W*d^2))));  
RHS = RHSsmallPsi + RHSbigPsi;

%As per Rusak, Wang 
%if the RHS is negative - set it to 0
RHS(RHS < 0) = 0;
wrstar = (Psi(2) - Psi(1))/dr;
res(2:end-2) = LHS - RHS;
res(1) = Psi(1);
res(end-1) = Psi(end) - 0.5*W*R^2;
res(end) = wrstar;
%%in case i want to see the full error
%dot(res,res)

end