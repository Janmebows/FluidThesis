%% Plot solutions to the steady axisymmetric flow in cylindrical coordinates


close all
clear all

syms r R k W

Amat = [r*besselj(1,k*r),r*bessely(1,k*r);...
        R*besselj(1,k*R),R*bessely(1,k*R)];
b = [-W*r^2 / 2 ; 0];
symx = Amat\b;
A = matlabFunction(symx(1));
B = matlabFunction(symx(2));



%rstar: basic vortex breakdown
%rhat: rankine body
R = 1;
W = 1;
rstar = 0.4* R;
kguess = 1;

deter = @(r,k,R)besselj(1,k.*r).*bessely(1,k.*R) - besselj(1,k.*R).*bessely(1,k.*r);
[r,k] = meshgrid(linspace(0.01,R),linspace(0.01,10));
Atimesdet =@(r,k,R) -0.5*W*r.*bessely(1,k*R) ;
Btimesdet =@(r,k,R) 0.5*W*r .*besselj(1,k*R);
wrstartimesdet = @(r,k,R) W*deter(r,k,R) + k.*(Atimesdet(r,k,R).*besselj(0,k.*r) + Btimesdet(r,k,R).*bessely(0,k.*r));
%wrstartimesdet = @(r,k,R) W*deter(r,k,R) +  (k./(r.*2)).* (Atimesdet(r,k,R).*(besselj(0,k.*r) - besselj(2,k.*r)) + Btimesdet(r,k,R).*(bessely(0,k.*r) - bessely(2,k.*r)));

contour(k*R,r/R,wrstartimesdet(r,k,R),[0,0],'k')
xlabel('kR')
ylabel('r^*/R')
legend('$$W(r^*) =0$$','interpreter','latex')

    
%solve wrstar = 0
%for various r,k combinations


r = linspace(0.01,R,1000);


%specify rstar and it gives me k

W = 1;
R = 1;

[k,fval,exitflag,output,jacobian] =fsolve(@(k)wrstar(k,rstar,W,R),kguess,optimoptions('fsolve','OptimalityTolerance', 1e-5)); 

psi = @(r) 0.5*W*r.^2 + r.*(A(R,W,k,rstar).*besselj(1,k.*r) + B(R,W,k,rstar).*bessely(1,k.*r));
w =@(r) W + k.*(A(R,W,k,rstar).*besselj(0,k.*r) + B(R,W,k,rstar).*bessely(0,k.*r));
legendentries = "$$r^*$$ = " + num2str(rstar) + ", k = " +num2str(k);
figure
hold on

%Enforce the construction: r < rstar gives w = 0
wplot = w(r);
wplot(r<rstar) = 0;
plot(r,wplot)
xlabel("r")
ylabel("w")
grid on
%just a line through w= zero
plot([0,r(end)],[0,0],'-k')
legend(legendentries,'interpreter','latex')


function out = wrstar(k,r,W,R)
deter = @(r,k,R)(besselj(1,k.*r).*bessely(1,k.*R) - besselj(1,k.*R).*bessely(1,k.*r));
Atimesdet =@(r,k,R) -0.5*W*r .*bessely(1,k*R) ;
Btimesdet =@(r,k,R) 0.5*W*r.*besselj(1,k*R);
out = W*deter(r,k,R) + k.*(Atimesdet(r,k,R).*besselj(0,k.*r) + Btimesdet(r,k,R).*bessely(0,k.*r));
end


