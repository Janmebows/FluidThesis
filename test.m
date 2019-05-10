%Obtain solution regions to the rankine body problem
%
%

close all
clear all
W = 1;
r0 = 1;
R = 1;
npts = 100;
[rstar,rhat,k] = meshgrid(linspace(0,R,npts),linspace(0,R,npts),linspace(0,20,npts));
C = 0.5 * W * (R.^2 - r0.^2)/(R.^2 - rhat.^2 ) ;
Ad = (rstar.*bessely(1,k.*rstar)).*(2*C -W) +k.*bessely(0,k.*rhat).*(0.5*W*rstar.^2);
Bd = -rstar.*besselj(1,k.*rstar).*(2*C-W)-k.*besselj(0,k.*rhat).*(0.5*W*rstar.^2);
deter = k.*rstar.*(besselj(0,k.*rhat).*bessely(1,k.*rstar) - bessely(0,k.*rhat).*besselj(1,k.*rstar));

wind = W*deter + k.*(Ad.*besselj(0,k.*rstar) + Bd.*bessely(0,k.*rstar));
psid = 0.5*W*deter.*(rhat.^2-r0^2) + rhat.*(Ad.*besselj(1,k.*rhat) + Bd.*bessely(1,k.*rhat));
%p = patch(isosurface(rstar,rhat,k,wind,0));
%p.FaceColor = 'red';
%p.EdgeColor = 'none';
isosurface(rstar,rhat,k,wind,0)
camlight
hold on
%isosurface(rstar,rhat,k,psid,0)

%isosurface(rstar,rhat,k,wind.^2 + psid.^2,0.01)
hold off
legend
view([30,60])
%so that the axes don't stretch weirdly
%axis vis3d
xlabel("rstar")
ylabel("rhat")
zlabel("k")
