syms r delta z pi Z omegaB
assume(r,'positive')
phiB = r * besselj(1,2*r*omegaB);

psi = 0.5*r^2 + delta * phiB * sin(pi*z/(2*Z));
lhs = simplify(diff(psi,z,2) + diff(psi,r,2) - diff(psi,r)/r);

phiByy = - 4 * omegaB^2 * phiB/ r^2;
eta = - delta * r *(phiByy - (pi/(2*Z))^2 * phiB/r^2) * sin(pi * z / (2*Z));
syms y omega phiB
v = (1./sqrt(2*y)).*((2*omega*y) + (2*delta*omega*phiB.*sin(pi*z/(2*Z))));
v = simplify( subs(v,y,r^2/2))
rhs = simplify(-r * eta);
"lhs"
pretty(lhs)


"rhs"
pretty(rhs)
