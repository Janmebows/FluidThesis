function [psiEq] = SymSquireLong1D(u,v,w,psi)
syms r
C = r*v;
dPsi = diff(psi,r);
dC = dPsi * diff(C,r);
dH = dPsi * (0.5*diff(u^2+v^2+w^2,r)+C^2/r^3);
rhs = simplify(r^2*dH - C*dC);
syms Psi(r)
psiEq = diff(Psi,r,2) - diff(Psi,r)/r ==  rhs;


end