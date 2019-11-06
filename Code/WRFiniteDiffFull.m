close all 
clear all
%parameters
nPtsR = 60;
nPtsZ = 100;
R = 1;
Z = 2.5;
W = 1;
Omega = 2;
dt = 0.05;
tEnd = 6;
%upstream flow
PsiInit = @(r) 0.5*W*r.^2;
VInit = @(r) Omega*r;
%put them in a params struct
params.nPtsR = nPtsR;
params.nPtsZ = nPtsZ;
params.R = R;
params.Z = Z;
params.W = W;
params.Omega = Omega;
r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
[rmat, zmat] = ndgrid(r,z);
dr = r(2)-r(1);
dz = z(2)-z(1);
%initial conditions
psi = PsiInit(rmat);
%using eta = -dwdr = - ddr(ddr(psi)./rmat)
eta = zeros(size(rmat)); %temporary
v   = VInit(rmat);


%more readable way for BCs
rIsZero = rmat==0;
zIsZero = zmat==0;
rIsR = rmat==R;
zIsZ = zmat==Z;

contours = linspace(0,0.1,20);
%preallocate matrices
etaStar = zeros(size(rmat));
vStar = zeros(size(rmat));
%precompute LU decomp for A * Psi = -r * eta
[L,U,A] = solveAAndDecompose(params);
figure
for t=0:dt:tEnd
    %steps 3,4
    dpsidz = ddz(psi,z);
    dpsidr = ddr(psi,r);
    J = @(x) dpsidz .* (ddr(x,r)) - (dpsidr .* (ddz(x,z)));
    %dvdt is exploding
    dvdt = J(v)./rmat + v.*dpsidz./rmat.^2 ;
    temp  = eta./rmat;
    %handle the 0/0 problem?
    temp(eta==0) =0;
    
    temp(1,:) = temp(2,:);
    detadt =  2*v.*ddz(v,z)./rmat;
    detadt(rIsZero)= 0;
    detadt = J(temp) +detadt;
    
        %step 7
    %iteration process
    %first iteration t=0 is odd step
%     if(mod(t/2+1,dt)==0)
%         etaStar = eta + dt .* detadt;
%         vStar = v + dt.* dvdt;
%         
%         dvStardt = J(vStar)./rmat + v.*dpsidz./rmat.^2 ;
%         detaStardt = J(etaStar./rmat) + 2*v.*ddz(vStar,z)./rmat;
%         eta = eta + dt.*detaStardt;
%         v = v + dt.*dvStardt;
%     else
          eta = eta + dt*detadt;
          v = v + dt*dvdt;
%     end

    %r=0
    v(rIsZero) = 0;
    %r=R
    v(rIsR) = VInit(R);
    %z=0
    v(zIsZero) = VInit(r);
    %z=Z
    v(zIsZ) = v(:,end-1);
    
    
    
        %step 5
    rhs = - rmat.*eta;
    
    %%PSI BCs----
    %r=0
    rhs(rIsZero) = 0;
    %r=R
    rhs(rIsR) = PsiInit(R);
    %z=0
    rhs(zIsZero) =PsiInit(r);
    %z=Z
    rhs(zIsZ) =rhs(:,end-1);
    

    
    
    rhsvec = rhs(:);
    y = L\rhsvec;
    psiV= U\y;
        %break if numerics break
    if any(isnan(psiV))
        "failed"
        break
    end
    
    psi = Vec2Mat(psiV,nPtsR,nPtsZ);
    

        
    %eta bcs 
    %r=0
    eta(rIsZero) = 0 ;
    %r=R
    %OLD
    %eta(end,:) = -2 * psi(end-1,:)/(R*dr^2);
    eta(rIsR) = (1/R).*((psi(end,:) - psi(end-1,:))/dr) - (1/R) .*((psi(end,:) - 2*psi(end-1,:) + psi(end-2,:))/(dr^2));
    %z=0
    %not sure how to handle this one...
    %OLD
    %eta(:,1) = -2* psi(:,2)./(r'*dz^2);
    %this line may be slightly off?
    i = 2:length(r)-1;
    eta(i,1) = (1./r(i).^2)'.*((psi(i+1,1) - psi(i-1,1))/dr) ...
        - (1./r(i))' .*((psi(i+1,1) - 2*psi(i,1) + psi(i-1,1))/(dr^2));
    %z=Z
    %OLD
    %eta(:,end) = -2*psi(:,end-1)./(r'*dz^2);
    eta(zIsZ) = eta(:,end-1);
    
    
    
    %display result
    surf(rmat,zmat,psi)
    axis([0,R,0,Z,0,inf])
    
    %contour(rmat,zmat,psi,contours)
    %contour(zmat,rmat,psi',50)
    %contourf(zmat,rmat,psi,20)
    %caxis([0,1])
    %colorbar
    xlabel("r")
    ylabel("z")
    zlabel("psi")
    title("t = "+ t)
    drawnow
    pause(0.05)

    


end




function [L,U,A] = solveAAndDecompose(params)
%Get the matrix A and decompose it into L U
%obtain the LU factorisation of the 
%finite difference formula for
%psi_zz + psi_rr - psi_r / r = A*psi
%tested and working
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

function ind = cvtIndex(i,j,nPtsR)
%converts an i,j into a single index ind 
%for finding psi_{i,j} stored in psiV(ind)
    ind = i+(j-1)*nPtsR;
end

%convert vector to matrix (tested)
function mat = Vec2Mat(vec,ndim,mdim)
    mat = zeros(ndim,mdim);
    for j=1:mdim
        lind = cvtIndex(1,j,ndim);
        rind = cvtIndex(0,j+1,ndim);
        mat(:,j) = vec(lind:rind);
    end
end

%central differences derivative wrt z (tested)
function dxdz = ddz(x,z)
    dxdz = zeros(size(x));
    dz = z(2) - z(1);
    j = 2:length(z)-1;
    dxdz(:,j) = (x(:,j+1) - x(:,j-1))/(2*dz);
end


%central differences derivative wrt r (tested)
function dxdr = ddr(x,r)

    dxdr = zeros(size(x));
    dr = r(2) - r(1);
    i = 2:length(r)-1;
    dxdr(i,:) = (x(i+1,:) - x(i-1,:))/(2*dr);
end