clear all 
close all
R = 1;
rdag = 0.9*R;
W = 1;
kguess = 0.8;


A = @(k,rdag) (-W*rdag)./(2*besselj(1,k.*rdag));
w = @(k,r,rdag) W + k.*A(k,rdag).*besselj(0,k.*r);


[r,k] = meshgrid(linspace(0.01,R));
contour(k*R,r/R,w(k,r,r),[0,0],'k')
xlabel('kR')
ylabel('$$r^\dagger/R$$','interpreter','latex')
legend('$$W(r^\dagger) =0$$','interpreter','latex')

k =fsolve(@(k) w(k,rdag,rdag),kguess,optimoptions('fsolve','OptimalityTolerance', 1e-8));

r = linspace(0,R);
figure
plot(r,w(k,r,rdag))
hold on
plot([0,R],[0,0])
plot([rdag,rdag],[w(k,r(1),rdag),w(k,r(end),rdag)])
%xlabel("r")
%ylabel("w(r)")
