syms u v w r g d C H psi R
u = 0
v = g/(2*pi*r) * (1-exp(-r^2/d^2))
C = r*v
rsub = sqrt(2*psi/w)
C = subs(C,r,rsub)
C = simplify(C)
dC = diff(C,psi)
dC = simplify(dC)

H = 1/2 * (u^2 + v^2 + w^2) + int(v^2/r,r)
H = subs(H,r,rsub)
H = simplify(H)
dH = diff(H,psi)
dH = simplify(dH)


dH = subs(dH,psi,0.5*w*r^2)
C = subs(C,psi,0.5*w*r^2)
dC = subs(dC,psi,0.5*w*r^2)
syms psi(r)
ode = diff(psi,r,2) - 1/r * diff(psi,r) == r^2 * dH - C* dC
dsolve(ode)
