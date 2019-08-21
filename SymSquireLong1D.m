function [psiEq,rhs] = SymSquireLong1D(u,v,w,psi)
syms r W 
C = r*v;
dPsi = diff(psi,r);
dC = dPsi * diff(C,r);
dH = dPsi * (0.5*diff(u^2+v^2+w^2,r)+C^2/r^3);
%psi = 0.5*W*r^2;
%this psi is different from before
syms Psi(r)
dC = subs(dC,r,sqrt(2*Psi/W))
dH = subs(dH,r,sqrt(2*Psi/W))
rhs = simplify(r^2*dH - C*dC);
%rhs = subs(rhs,r,sqrt(2*Psi/W))
%syms Psi(r)
psiEq = diff(Psi,r,2) - diff(Psi,r)/r ==  rhs;

end