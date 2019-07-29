syms u v w r g d C H Psi R W
u = 0;
v = g/(2*pi*r) * (1-exp(-r^2/d^2))
psi = 0.5*W*r^2;
w = W
[psiEq,rhs] = SymSquireLong1D(u,v,w,psi)

dsolve(psiEq)
