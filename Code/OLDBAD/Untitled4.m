W = 1;
G = 4;
d = 0.25;
R = 1;
rstar = 0.05;
points = 100;


i = 2:points-1;
numSteps = 100000;

%relaxation method
Psi = zeros(1,points);
Psi(1) = 0;
Psi(end) = 0.5*W*R^2;
r = linspace(0,R,points);
dr = r(2)-r(1);
tol = dr;

Psi(i) = 0.5*W*i.^2;
%f is the RHS of the SL equation
f = @(r,Psi) G^2/(2*W*d^2*pi^2) * ((r.^2*W./(2*Psi) -1).*(exp(-2*Psi/(W*d^2)) - exp(-4*Psi/(W*d^2))));
for currentStep =1:numSteps
Psi(i) = 0.2*(Psi(i+1).*(2-dr./r(i)) +Psi(i)+ Psi(i-1).*(2 + dr./r(i)) - (f(r(i),Psi(i))*(2*dr^2)));
%BCs
Psi(r-rstar < tol) = 0;
Psi(end) = 0.5*W*R^2;
end
w = (Psi(i+1) - Psi(i-1))./r(i);
subplot(2,1,1)
plot(r(i),w)
subplot(2,1,2)
plot(r,Psi)

axis([0,R,0,inf])