syms u v w r G d C H Psi(r) R W pi
u = 0;
v = G/(2*pi*r) * (1-exp(-r^2/d^2))
psi = 0.5*W*r^2;
w = W

C = r*v;
H = 0.5 * (u^2 + v^2 + w^2) + int(C^2/r^3,r) 
%dPsi = diff(psi,r);
%dC = diff(C,psi);

%dH = dPsi * (0.5*diff(u^2+v^2+w^2,r)+C^2/r^3);

%psi = 0.5*W*r^2;
%this psi is different from before
syms Psi
C =subs(C,r,sqrt(2*Psi/W))
H = subs(H,r,sqrt(2*Psi/W))
dC = diff(C,Psi)
dH = diff(H,Psi)
% dC = subs(dC,r,sqrt(2*Psi/W))
% dH = subs(dH,r,sqrt(2*Psi/W))
rhs = simplify(r^2*dH - C*dC)
%rhs = subs(rhs,r,sqrt(2*Psi/W))

limit(rhs,Psi,0,'right')
taylor(rhs,Psi,0,'Order',1)
func1 = matlabFunction(rhs)
func2 = matlabFunction(taylor(rhs,Psi,0,'Order',2))

%psiEq = diff(Psi,r,2) - diff(Psi,r)/r ==  rhs;


