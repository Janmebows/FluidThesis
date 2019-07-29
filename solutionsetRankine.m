%plot rstar rhat solution space
W = 1;
r0 = 1;
R = 1;
npts = 100;
det =@(rhat,rstar,k) rhat.*(bessely(1,k.*rhat).*besselj(1,k.*rstar) - besselj(1,k.*rhat).*bessely(1,k.*rstar));
A =@(rhat,rstar,k) 1/det(rhat,rstar,k) .* (rhat.*bessely(1,k.*rhat).*(-0.5.*W.*rstar) -bessely(1,k.*rstar).*(0.5.*W.*(r0.^2-rhat.^2)));
B =@(rhat,rstar,k) 1/det(rhat,rstar,k) .* (-rhat.*besselj(1,k.*rhat).*(-0.5.*W.*rstar) +besselj(1,k.*rstar).*(0.5.*W.*(r0.^2-rhat.^2)));
C =@(rhat,rstar,k) 0.5.*(W + k.*(A(rhat,rstar,k).*besselj(0,k.*rhat) + B(rhat,rstar,k).*bessely(0,k.*rhat)));
D =@(rhat,rstar,k) 0.5.*W.*R.^2 -C(rhat,rstar,k);
rhateqn =@(rhat,rstar,k) (A(rhat,rstar,k).*besselj(0,k.*rhat) + B(rhat,rstar,k).*bessely(0,k.*rhat)).*(rhat.^2 -1) - W/k .* (r0.^2 - R.^2 -1);
rstareqn =@(rhat,rstar,k) W + k.*(A(rhat,rstar,k).*besselj(0,k.*rstar) + B(rhat,rstar,k).*bessely(0,k.*rstar));


[rstar,rhat,k] = meshgrid(linspace(0,R,npts),linspace(0,R,npts),linspace(0,20,npts));

%p = patch(isosurface(rstar,rhat,k,wind,0));
%p.FaceColor = 'red';
%p.EdgeColor = 'none';
%isosurface(rstar,rhat,k,rhateqn(rhat,rstar,k),0)
camlight
hold on
%isosurface(rstar,rhat,k,rstareqn(rhat,rstar,k),0)

isosurface(rstar,rhat,k,rhateqn(rhat,rstar,k).^2 + rstareqn(rhat,rstar,k).^2,0.1)
hold off
legend
view([30,60])
%so that the axes don't stretch weirdly
%axis vis3d
xlabel("rstar")
ylabel("rhat")
zlabel("k")
