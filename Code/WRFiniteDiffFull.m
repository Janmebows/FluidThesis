close all 
clear all
nPtsR = 100;
nPtsZ = 100;
R = 1;
Z = 1.5;
W = 1;
params.nPtsR = nPtsR;
params.nPtsZ = nPtsZ;
params.R = R;
params.Z = Z;
params.W = W;
dt = 0.05;
tEnd = 10;
r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
[rmat, zmat] = ndgrid(r,z);
dr = r(2)-r(1);
dz = z(2)-z(1);
%initial conditions
psi = 0.5*W*rmat.^2;
eta = 0*rmat; 
v = W*zmat;
%preallocate matrices
etaStar = zeros(size(rmat));
vStar = zeros(size(rmat));
%precompute LU decomp for A * Psi = r * eta
[L,U,A] = solveAAndDecompose(params);
figure
for t=0:dt:5
    
    %steps 3,4
    dpsidz = ddz(psi,z);
    dpsidr = ddr(psi,r);
    J = @(x) dpsidz .* (ddz(x,z)) - dpsidr .* (ddr(x,r));
    dvdt = J(v)./rmat + v.*dpsidz./rmat.^2 ;
    temp  = eta./rmat;
    %handle the 0/0 problem
    temp(eta==0) =0;
    detadt = J(temp) + 2*v.*ddz(v,z)./rmat;
%     detadt(isnan(detadt)) = 0;
    %step 5
    rhs = - rmat.*eta;
    rhs = PsiBCs(rhs,params);
    rhsvec = rhs(:);
    y = L\rhsvec;
    psiV= U\y;
    psi = Vec2Mat(psiV,nPtsR,nPtsZ);

    %step 7
    %iteration process
    %first iteration t=0 is odd step
    if(mod(t/2+1,dt)==0)
        etaStar = eta + dt .* detadt;
        vStar = v + dt.* dvdt;
        
        dvStardt = J(vStar)./rmat + v.*dpsidz./rmat.^2 ;
        detaStardt = J(etaStar./rmat) + 2*v.*ddz(vStar,z)./rmat;
        eta = eta + dt.*etaStar;
        v = v + dt.*vStar;
    else
        eta = eta + dt.*eta;
        v = v + dt.*v;
    end
    
    %display result
    
    surf(rmat,zmat,psi)
    axis([0,R,0,Z,0,inf])
    
    %contourf(zmat,rmat,psi,20)
    %caxis([0,20])
    %colorbar
    xlabel("r")
    ylabel("z")
    zlabel("psi")
    drawnow
    pause(0.05)
    

    v(1,:) = 0;
    v(end,:) = 0;
    v(:,1) = r;
    v(:,end) = 0;
        
    %eta bcs 

    eta(end,:) = -2 * psi(end-1,:)/(dr^2);
    
    
    eta(:,1) = -2* psi(:,2)./(r'*dz^2);
    eta(:,end) = -2*psi(:,end-1)./(r'*dz^2);
    eta(1,:) = 0 ;



end




function [L,U,A] = solveAAndDecompose(params)
%obtain the LU factorisation of the 
%finite difference formula for
%psi_zz + psi_rr - psi_r / r = A*psi
R = params.R;
Z = params.Z;
nPtsR = params.nPtsR;
nPtsZ = params.nPtsZ;
r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
dr = r(2) - r(1);
dz = z(2) - z(1);
invdr = 1/dr;
invdz = 1/dz;
A = spalloc(nPtsR*nPtsZ,nPtsR*nPtsZ,nPtsR*nPtsZ + 5*(nPtsR-1)*(nPtsZ-1));

%loop over j and i
%i can't seem to vectorise this nicely so double loop for now
for j=2:nPtsZ-1
   for i=2:nPtsR-1
       %b = a + 1
       A(i+(j-1)*nPtsR,i+1+(j-1)*nPtsR) = invdr^2 - (invdr./(2*r(i)));
       %b = a
       A(i+(j-1)*nPtsR,i+(j-1)*nPtsR) = -2*(invdz^2 + invdr^2);
       %b = a - 1
       A(i+(j-1)*nPtsR,i-1+(j-1)*nPtsR) = invdr^2 + (invdr./(2*r(i)));
       %b = a - m
       A(i+(j-1)*nPtsR,i+(j-2)*nPtsR) = invdz^2;
       %b = a +m
       A(i+(j-1)*nPtsR,i+j*nPtsR) = invdz^2; 
   end    
end

%BCs
for i=1:nPtsR
    %z = 0 => j=1
    ind = cvtIndex(i,1,nPtsR);
    A(ind,ind) = 1;
    %z = Z => j = nPts
    ind = cvtIndex(i,nPtsZ,nPtsR);
    A(ind,ind) = 1;

end   
for j=1:nPtsZ
    ind = cvtIndex(1,j,nPtsR);
    A(ind,ind) = 1;
    %r=R
    ind = cvtIndex(nPtsR,j,nPtsR);
    A(ind,ind) = 1;
end   
[L, U] = lu(A);
end
% plot3(r,z,eta)
function ind = cvtIndex(i,j,nPtsR)
    ind = i+(j-1)*nPtsR;
end

function mat = Vec2Mat(vec,ndim,mdim)
    mat = zeros(ndim,mdim);
    for j=1:mdim
        lind = cvtIndex(1,j,ndim);
        rind = cvtIndex(0,j+1,ndim);
        mat(:,j) = vec(lind:rind);
    end
end
function rhs = PsiBCs(rhs,params)
%     W = params.W;
%     R = params.R;
    %%BCs----
    %r=0
    rhs(1,:) = 0;
    %r=R
    rhs(end,:) = 0;
    %z=0
    rhs(:,1) =0;
    %z=Z
    rhs(:,end) =0;

end


function dxdz = ddz(x,z)
    dxdz = zeros(size(x));
    dz = z(2) - z(1);
    j = 2:length(z)-1;
    dxdz(:,2:end-1) = (x(:,j+1) - x(:,j-1))/(2*dz);
end

function dxdr = ddr(x,r)

    dxdr = zeros(size(x));
    dr = r(2) - r(1);
    i = 2:length(r)-1;
    dxdr(2:end-1,:) = (x(i+1,:) - x(i-1,:))/(2*dr);
end