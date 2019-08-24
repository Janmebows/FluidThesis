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



func = @(eta,Psi,rstar)psiODE(eta,Psi,rstar,R,W,k);

step = 0.1;
errtol = 1e-10;
lastDir = 0;
flag = false;
%I am making a big assumption on the relationship between the
% size of rstar and the error
while ~flag
    [r,Psi] = ode45(@(r,Psi)func(r,Psi,rstarguess),[rstarguess,R],[0,0]);
    
    if rstarguess < 0 || rstarguess > 1
        fprintf("Failed\n\n")
        return
    end
    err = Psi(end,1) - 0.5*W*R^2;
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
end

rstar = rstarguess

w = Psi(:,2)./r;
plot(r,Psi(:,1))
ylabel("$$\Psi$$")
axis([0,R,0,inf])
figure
plot(r,w)
ylabel("$$w$$")
axis([0,R,0,inf])





function out= psiODE(r,Psi,rstar,R,W,k)
rhs = (k^2/(2*W)).* r.^2 - k^2 .* Psi(1,:); 
dPsidr = Psi(2,:);
d2Psidr2 =dPsidr./r + rhs;
out = [dPsidr ; d2Psidr2];

end