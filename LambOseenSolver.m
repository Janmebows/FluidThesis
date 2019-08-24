close all
R = 1;
W = 1;
G = 4.5;
d = 0.5;
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
func = @(eta,Psi,rstar)psiODELambOseen(eta,Psi,rstar,R,W,G,d);

%flag is false when a solution hasn't been found
flag = 0;

%solution stepping
step = 0.1;
errtol = 1e-10;
lastDir = 0;

%I am making a big assumption on the relationship between the
% size of rstar and the error
rstarguess=  0.5;
step = 0.1;
lastDir = 0;

func = @(eta,Psi,rstar)psiODELambOseen(eta,Psi,rstar,R,W,G,d);
 
flag = 0;
while flag==0
    [r,psi] = ode45(@(r,Psi)func(r,Psi,rstarguess),[rstarguess,R],[0,0]);
    
    if rstarguess < 0 || rstarguess > 1
        %fprintf("Failed\n\n")
        flag = -1;
    end
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
        flag = 1;
    end 
end

if flag ==1
rstar = rstarguess
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


function out= psiODELambOseen(r,Psi,rstar,R,W,G,d)

GamTerm = G^2/(2*W*d^2*pi^2);
expTerm = -2/(W*d^2);
  if Psi(1,:) < 0.0001
      PsiSeries = 1 - (-1+Psi(1,:)) + (-1+Psi(1,:)).^2 -(-1+Psi(1,:)).^3;
      %rhs = GamTerm* ((r.^2*W .*PsiSeries/2) -1).*(exp(Psi(1,:)*expTerm ) - exp(2*Psi(1,:)*expTerm )); 
      rhs = (G^2*r^2)/(2*W*d^4*pi^2);
  else
     rhs = GamTerm* ((r.^2*W ./(2*Psi(1,:))) -1).*(exp(Psi(1,:)*expTerm ) - exp(2*Psi(1,:)*expTerm ));  
  end
dPsidr = Psi(2,:);
d2Psidr2 = dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end