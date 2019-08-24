close all
R = 1;
W = 1;
G = 10;
d = 0.8331;
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


%make it easier to call
func = @(eta,Psi,rstar)psiODERotating(eta,Psi,rstar,R,W,G,d);

%flag is false when a solution hasn't been found
flag = false;

%solution stepping
step = 0.1;
errtol = 1e-10;
lastDir = 0;

%I am making the fatal assumption that overshooting implies a 
while ~flag
    [r,psi] = ode45(@(r,Psi)func(r,Psi,rstarguess),[rstarguess,R],[0,0]);
    err = psi(end,1) - 0.5*W*R^2;
    if err > errtol
        fprintf("r^* is too small \n")
        rstarguess = rstarguess + step;
        %decrease step size if we change direction
        if lastDir ==1
            continue;
        elseif lastDir ==-1
            step = step/10;
        else 
            lastDir =1;
        end
    elseif err < -errtol
        fprintf("r^* is too big \n")
        rstarguess = rstarguess - step;
        if lastDir ==-1
            continue;
        elseif lastDir ==1
            step = step/10;
        else 
            lastDir =-1;
        end
    else 
        flag = true
    end 
    if rstarguess < 0 || rstarguess > 1
        fprintf("Failed\n\n")
        return
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

function out= psiODERotating(r,Psi,rstar,R,W,G,d)
k = G;
rhs = (k^2/(2*W)).* r.^2 - k^2 .* Psi(1,:); 
dPsidr = Psi(2,:);
d2Psidr2 =dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end
function out= psiODELambOseen(r,Psi,rstar,R,W,G,d)
% rhs = (k^2/(2*W)).* r.^2 - k^2 .* Psi(1,:); 
% dPsidr = Psi(2,:);
% d2Psidr2 =dPsidr./r + rhs;
% out = [dPsidr ; d2Psidr2];
GamTerm = G^2/(2*W*d^2*pi^2);
expTerm = -2/(W*d^2);
  if Psi(1,:) < 0.0001
%      PsiSeries = 1 - (-1+Psi(1,:)) + (-1+Psi(1,:)).^2 -(-1+Psi(1,:)).^3;
%      rhs = GamTerm* ((r.^2*W .*PsiSeries) -1).*(exp(Psi(1,:)*expTerm ) - exp(2*Psi(1,:)*expTerm )); 
    rhs = (G^2*r^2)/(2*W*d^4*pi^2);
  else
     rhs = GamTerm* ((r.^2*W ./(2*Psi(1,:))) -1).*(exp(Psi(1,:)*expTerm ) - exp(2*Psi(1,:)*expTerm ));  
  end
dPsidr = Psi(2,:);
d2Psidr2 = dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end