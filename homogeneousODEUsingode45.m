close all
R = 1;
W = 1;
k = 5;
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


%singular de must include the singular point
%in this case dpsi/r(0) = 0
initguess = @(eta) [0.5*W*eta^2,W*eta]%[1.5*eta.^2+1,1.5*eta.^2];
%may need to find a way to solve right to left rather than left to right.
numPts = 100;
solinit=bvpinit(linspace(0,1,numPts),initguess,rstarguess);
%options =bvpset('RelTol',1e-5,'AbsTol',1e-8,'Nmax',1e6,'Vectorized','on');%,'SingularTerm',[0,0;0,1]);
options =bvpset('RelTol',1e-5,'AbsTol',1e-8,'Nmax',1e6);%,'Vectorized','on');%,'SingularTerm',[0,0;0,1]);
func = @(eta,Psi,rstar)psiODE(eta,Psi,rstar,R,W,k);
%bound = @(psirstar,psiR,rstar)boundaries(psirstar,psiR,rstar,R,W);
%sol=bvp5c(func,bound,solinit,options);
%if too big
flag = false;
step = 0.00005;
errtol = 1e-5;
while ~flag
    [r,psi] = ode45(@(r,Psi)func(r,Psi,rstarguess),[rstarguess,R],[0,0]);
    err = psi(end,1) - 0.5*W*R^2
    if err > errtol
    rstarguess = rstarguess + step
    elseif err < -errtol
        rstarguess = rstarguess - step
    else 
        flag = true
    end 
    if rstarguess < 0 || rstarguess > 1
        fprintf("Won't work\n")
        break;
    end
    
end
rstar = rstarguess

w = psi(:,2)./r;
plot(r,psi(:,1))
ylabel("$$\Psi$$")
axis([0,R,0,inf])
figure
plot(r,w)
ylabel("$$w$$")
axis([0,R,0,inf])



function res = boundaries(psirstar,psiR,rstar,R,W)
%the boundary conditions
%psi(r*) = 0
%psi(R) = 0.5*W*R^2
%w(r*) = dPsidr/r(r*) =0
res = [psirstar(1);
       psiR(1) - 0.5*W*R^2;
        psirstar(2)/rstar];
end


function r = etaTor(eta,rstar,R)
%eta is 0,...,1
%convert to r = rstar, ..., R
r = eta.*(R-rstar) + rstar;

end

function eta = rToeta(r,rstar,R)
eta = (r-rstar)./(R-rstar);
end

function out= psiODE(r,Psi,rstar,R,W,k)
rhs = (k^2/(2*W)).* r.^2 - k^2 .* Psi(1,:); 
dPsidr = Psi(2,:);
d2Psidr2 =dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end