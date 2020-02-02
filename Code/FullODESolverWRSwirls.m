close all 
clear all
%parameters
nPtsR = 20;
nPtsZ = 50;
R = 1;
Z = 6;
W = 1;
Omega = 1.9; %is this the same omega as ours?
%omegaB is the first zero of J_1 which is near 3
omegaB = fsolve(@(x) besselj(1,x),3)/2; %good
tEnd = 6;
delta = -1.111;
n=1; %i think they use n=1 throughout without explicitly saying it
     %so ill use that
%upstream flow
%something wrong here it should be 0 when y = 1/2
phiB = @(y) sqrt(2*y).*besselj(1,2*omegaB.*sqrt(2*y)); %good
omegamn = sqrt(omegaB^2 + (2*n-1)^2 * pi^2/(16*Z^2));
phiByy = @(y) -(4*omegamn^2 - ((2*n-1)^2 *pi^2)/(4*Z^2))*phiB(y)./(2*y);%TEMP
PsiInit = @(r,z) r.^2/2 + delta *phiB(r.^2/2).*sin(pi*z/(2*Z));
etaInit = @(r,z) r.*(-delta*(phiByy(r.^2/2) - (pi/(2*Z))^2.*phiB(r.^2/2)./(2*r.^2/2))).*sin(pi*z/(2*Z));
VInit = @(r,z) (1./r).*(2*Omega*r.^2/2 + 2*delta*Omega*phiB(r.^2/2).*sin(pi*z/(2*Z)));
%put them in a params struct
params.nPtsR = nPtsR;
params.nPtsZ = nPtsZ;
params.R = R;
params.Z = Z;
params.W = W;
params.Omega = Omega;
params.PsiInit = PsiInit;
params.VInit = VInit;
params.etaInit = etaInit;
params.psi = [];
r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
params.r = r;
params.z = z;
[rmat, zmat] = ndgrid(r,z);
dr = r(2)-r(1);
dz = z(2)-z(1);
%initial conditions
psi = PsiInit(rmat,zmat);
%using eta = -dwdr = - ddr(ddrr(psi)./rmat)
eta = etaInit(rmat,zmat);
v   = VInit(rmat,zmat);



contours = linspace(0.01,0.49,25);
%preallocate matrices
%etaStar = zeros(size(rmat));
%vStar = zeros(size(rmat));
%precompute LU decomp for A * Psi = -r * eta
[L,U,A] = solveAAndDecompose(params);
params.L = L;
params.U = U;

etaVInit(:,:,1) = eta(2:end-1,2:end-1);
etaVInit(:,:,2) = v(2:end-1,2:end-1);

t = linspace(0,6,100);
[t,out] = ode45(@(t,in)DE(t,in,params),t,etaVInit);

%post process the data

eta = zeros(nPtsR,nPtsZ,length(t));
v = zeros(nPtsR,nPtsZ,length(t));

%unfortunately i've had to do this in a 
%loop since reshape isn't cooperating
for i = 1:length(t)
    temp = reshape(out(i,:), [nPtsR-2,nPtsZ-2,2]);
    eta(2:end-1,2:end-1,i) = temp(:,:,1);
    v(2:end-1,2:end-1,i) = temp(:,:,2);
    
        
    %psi
    rhs = - rmat.*eta(:,:,i);
    
    %%PSI BCs----
    %r=0
    rhs(1,:) = 0;
    %r=R
    rhs(end,:) = PsiInit(R,z);
    %z=0
    rhs(:,1) =PsiInit(r,0);
    %z=Z -> z derivative is zero
    rhs(:,end) =0;
    

    
    
    rhsvec = rhs(:);
    y = L\rhsvec;
    psiV= U\y;        %break if numerics break
    if any(isnan(psiV))
        error("failed")
    end
    
    psi(:,:,i) = Vec2Mat(psiV,nPtsR,nPtsZ);
    
    
    eta(:,:,i) = EtaBCs(eta(:,:,i),psi(:,:,i),params);
    v(:,:,i) = VBCs(v(:,:,i),params);

end
 
%%plotting 
figure

for i = 1:length(t)
    contour(rmat,zmat,psi(:,:,i),contours)
    %contour(zmat,rmat,psi',50)
    %contourf(zmat,rmat,psi,20)
    %caxis([0,1])
    colorbar
    xlabel("r")
    ylabel("z")
    zlabel("psi")
    title("t = "+ t(i))
    drawnow
    %Frame(i) = getframe(gcf);
    pause

    
end
%%
%%%save to file
 vwriter = VideoWriter('contourseta.mp4','MPEG-4');
 open(vwriter)
 writeVideo(vwriter,Frame);
 close(vwriter)
function out = DE(~,in,params)
%%%%%unpack variables

R = params.R;
Z = params.Z;
nPtsR = params.nPtsR;
nPtsZ = params.nPtsZ;
%format in vector into nPtsR-2,nPtsZ-2,2 matrix
in = reshape(in,[nPtsR-2,nPtsZ-2,2]);
eta = zeros(nPtsR,nPtsZ);
v = zeros(nPtsR,nPtsZ);
eta(2:end-1,2:end-1) = in(:,:,1);
v(2:end-1,2:end-1) = in(:,:,2);

PsiInit = params.PsiInit;
VInit = params.VInit;
L = params.L;
U = params.U;

%might pass these in as params
r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
[rmat, zmat] = ndgrid(r,z);

dr = r(2)-r(1);
% dz = z(2)-z(1);

%more readable way for BCs
rIsZero = rmat==0;
zIsZero = zmat==0;
rIsR = rmat==R;
zIsZ = zmat==Z;
%%%%%
%solve for psi
    rhs = - rmat.*eta;
    
    %%PSI BCs----
    %r=0
    rhs(1,:) = 0;
    %r=R
    rhs(end,:) = PsiInit(R,z);
    %z=0
    rhs(:,1) =PsiInit(r,0);
    %z=Z -> z derivative is zero
    rhs(:,end) =0;
    

    
    
    rhsvec = rhs(:);
    y = L\rhsvec;
    psiV= U\y;        %break if numerics break
    if any(isnan(psiV))
        error("failed")
    end
    
    psi = Vec2Mat(psiV,nPtsR,nPtsZ);
    

   
%     %r=0
%     v(rIsZero) = 0;
%     %r=R
%     v(rIsR) = VInit(R);
%     %z=0
%     v(zIsZero) = VInit(r);
%     %z=Z
%     v(zIsZ) = v(:,end-1);
%     
    v = VBCs(v,params);
    
        %step 5

        eta = EtaBCs(eta,psi,params);

    
     
%steps 3,4
    dpsidz = ddz(psi,z);
    dpsidr = ddr(psi,r);
    J = @(x) (dpsidz .* (ddr(x,r))) - (dpsidr .* (ddz(x,z)));
    %dvdt is exploding
    dvdt = J(v)./rmat + v.*dpsidz./rmat.^2 ;
    temp = eta./rmat;
    %handle the 0/0 problem?
    temp(eta==0) =0;
    
    temp(1,:) = temp(2,:);
    detadt =  2*v.*ddz(v,z)./rmat;
    detadt(rIsZero)= 0;
    detadt = J(temp) +detadt;
    
    out(:,:,1) = detadt(2:end-1,2:end-1);
    out(:,:,2) = dvdt(2:end-1,2:end-1);

    out = out(:);


end


function eta = EtaBCs(eta,psi,params)
r = params.r;
R = params.R;
z = params.z;
etaInit = params.etaInit;

dr = r(2)-r(1);
    %r=0
    eta(1,:) = 0;
    %r=R
    %this guy is problematic
    eta(end,:) = (1/R^2).*((psi(end,:) - psi(end-1,:))/dr) - (1/R) .*((psi(end,:) - 2*psi(end-1,:) + psi(end-2,:))/(dr^2));

    %z=0
    %not enforcing anything
    i = 2:length(r)-1;
    eta(i,1) = ((1./r(i).^2)'.*((psi(i+1,1) - psi(i-1,1))/(2*dr))) ...
        - ((1./r(i))' .*((psi(i+1,1) - 2*psi(i,1) + psi(i-1,1))/(dr^2)));
    %z=Z
    %OLD
    eta(:,end) = eta(:,end-1);
end

function v = VBCs(v,params)
r = params.r;
z = params.z;
R = params.R;
VInit = params.VInit;
    %r=0
    v(1,:) = 0;
    %r=R
    v(end,:) = VInit(R,z);
    %z=0
    v(:,1) = VInit(r,0);
    %z=Z
    v(:,end) = v(:,end-1);
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
    ind2 = cvtIndex(i,nPtsZ-1,nPtsR);
    A(ind,ind2) = -1;

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