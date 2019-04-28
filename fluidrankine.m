clear all
close all

%can choose
r0 = 0.2;
%rstar = 0.3;
%rhat = 0.5;
R = 1;
k = 1;
W = 1;


%rhat = sqrt(r0^2*(R^2 + 1)/(R^2 - 1  + 2*r0^2));
%C = 
%D = 
%A = @(rstar,rhat,k)
%w = @(rstar,rhat,k) W + k.*(A.*besselj(0,k.*rstar) + B.*bessely(0,k.*rstar));
%need to fix this
%its still the old code

C =@(rhat) 0.5 * W * (R^2 - r0^2)/(R^2 - rhat^2 );
D =@(rhat) 0.5 * W * R^2 * (r0^2 - rhat^2)/(R^2 - rhat^2);
psiouter = @(r,rhat) C(rhat)*r.^2 + D(rhat);

rhat = fsolve(@(r) psiouter(r,r) - 0.5*W*r0^2,0.2);

r = linspace(0,R);  
plot(r,psiouter(r,rhat))
hold on
plot([0,R],[0.5*W*r0^2,0.5*W*r0^2])


syms rstar

Amat = [k*besselj(0,k*rhat),k*bessely(0,k*rhat);...
        rstar*besselj(1,k*rstar),rstar*bessely(1,k*rstar)];
b = [2*C(rhat) - W;-W*rstar^2 / 2];
x = Amat\b;
A = x(1);
B = x(2);

%now obtain rhat
%using w(r_*) = 0
A = matlabFunction(A);
B = matlabFunction(B);
wout =2*C(rhat);

%    [rstar,fval,exitflag,output,jacobian] =fsolve(wout,0.5);

if rhat < rstar || exitflag <1
    win = @(r) W + k.*(A(rhat).*besselj(0,k.*r) + B(rhat).*bessely(0,k.*r));
    [rhat,fval,exitflag,output,jacobian] =fsolve(win,0.1);
    if exitflag < 1
        error("Doesn't work")
    end

end

