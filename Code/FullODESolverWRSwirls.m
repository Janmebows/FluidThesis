close all 
clear all
%%%%%%%%%%%
%parameters
%%%%%%%%%%%

omega = 1.9;
%omega = 1.95;
delta = -1.111;
%delta = 2.5;


nPtsR = 40;
nPtsZ = 100;
nPtst = 50;
R = 1;
Z = 6;
tEnd = Z^3 * 2;

%upstream flow

%omegaB is the first zero of J_1 which is near 3
omegaB = fsolve(@(x) besselj(1,x),3)/2;



phiB = @(r) r .* besselj(1,2*r*omegaB);
phiByy =@(r) - 4 * omegaB^2 * phiB(r)./ r.^2;
PsiInit = @(r,z) 0.5*r.^2 + delta * phiB(r) .* sin(pi*z/(2*Z));
EtaInit = @(r,z) - delta * r .*(phiByy(r) - (pi/(2*Z))^2 * phiB(r)./r.^2) .* sin(pi * z / (2*Z));
VInit = @(r,z) omega./r .*(r.^2 + 2 * delta * phiB(r) .* sin(pi*z/(2*Z)));



%%%verify that solid body works
%%%I.e. no perturbation
% PsiInit = @(r,z) 0.5 * r.^2;
% VInit = @(r,z) omega * r;
% etaInit = @(r,z) 0;


r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
[rmat, zmat] = ndgrid(r,z);


%put them in a params struct
params.nPtsR = nPtsR;
params.nPtsZ = nPtsZ;
params.R = R;
params.Z = Z;
params.omega = omega;
params.PsiInit = PsiInit;
params.VInit = VInit;
params.EtaInit = EtaInit;
params.psi = [];
params.r = r;
params.z = z;
params.rmat = rmat;
params.zmat = zmat;
dr = r(2)-r(1);
dz = z(2)-z(1);


%apply initial conditions
psi = PsiInit(rmat,zmat);
eta = EtaInit(rmat,zmat);
v   = VInit(rmat,zmat);


%precompute LU decomp for A * Psi = -r * eta
[L,U,A] = solveAAndDecompose(params);
params.L = L;
params.U = U;
params.psi = psi;

etaVInit(:,:,1) = eta(2:end-1,2:end-1);
etaVInit(:,:,2) = v(2:end-1,2:end-1);

t = linspace(0,tEnd,nPtst*tEnd);
[t,out] = ode45(@(t,in)DE(t,in,params),t,etaVInit);

%post process the data

eta = zeros(nPtsR,nPtsZ,length(t));
v = zeros(nPtsR,nPtsZ,length(t));

%unfortunately i've had to do this in a 
%loop since reshape isn't cooperating
for i = 1:length(t)-1
    temp = reshape(out(i,:), [nPtsR-2,nPtsZ-2,2]);
    eta(2:end-1,2:end-1,i) = temp(:,:,1);
    v(2:end-1,2:end-1,i) = temp(:,:,2);
        
    eta(:,:,i) = EtaBCs(eta(:,:,i),psi(:,:,i),params);
    v(:,:,i) = VBCs(v(:,:,i),params);
    psi(:,:,i+1) = PsiFromEta(eta(:,:,i), params);

end


%%plotting 
%%
psicontours = linspace(0.01,0.49,25);

for i = 1:floor(length(t)/100):length(t)
    fig1 = figure(1);
    contour(zmat',rmat',psi(:,:,i)',psicontours)
    colorbar
    xlabel("z")
    ylabel("r")
    title("\psi(t), t = "+ t(i))
    drawnow
    %Frame1(i) = getframe(figure(1));
    fig2=figure(2);
    contour(zmat',rmat',eta(:,:,i)')
    colorbar
    xlabel("z")
    ylabel("r")
    title("\eta(t), t = "+ t(i))
    drawnow
    
    
    %Frame2(i) = getframe(figure(2));
    
    
    fig3 =figure(3);
    contour(zmat',rmat',v(:,:,i)')
    colorbar
    xlabel("z")
    ylabel("r")
    zlabel("v")
    title("v(t), t = "+ t(i))
    drawnow
    %Frame3(i) = getframe(figure(3));
    
    %pause(0.0001)

    
end
%%
%%%save to file
% names = ["contourspsi.mp4","contourseta.mp4","contoursv.mp4"];
% 
% vwriter = VideoWriter("contourspsi.mp4",'MPEG-4');
% open(vwriter)
% writeVideo(vwriter,Frame1);
% close(vwriter)
% 
% vwriter = VideoWriter("contourseta.mp4",'MPEG-4');
% open(vwriter)
% writeVideo(vwriter,Frame2);
% close(vwriter)
% 
% vwriter = VideoWriter("contoursv.mp4",'MPEG-4');
% open(vwriter)
% writeVideo(vwriter,Frame3);
% close(vwriter)
% 
%   
%   
  
  
 function psi = PsiFromEta(eta, params)

rmat = params.rmat;
R = params.R;
z = params.z;
r = params.r;
L = params.L;
U = params.U;
nPtsR = params.nPtsR;
nPtsZ = params.nPtsZ;
PsiInit = params.PsiInit;
 %psi
    rhs = - rmat.*eta;
    
    %%PSI BCs----
    %r=0
    rhs(1,:)    = 0;
    %r=R
    rhs(end,:)  = PsiInit(R,z);
    %z=0
    rhs(:,1)    = PsiInit(r,0);
    %z=Z -> z derivative is zero
    rhs(:,end)  = 0;
    

    
    
    rhsvec = rhs(:);
    y = L\rhsvec;
    psiV= U\y;    

    psi = Vec2Mat(psiV,nPtsR,nPtsZ);
 
 end
function out = DE(~,in,params)
%%%Process:
%1. get psi everywhere for time t
%2. get inner points for eta and v for time t+deltat
%3. get psi everywhere for time t+deltat
%4 get BCs for eta and v for time t+deltat 
%out detadt and dvdt


%%%%%unpack variables
nPtsR = params.nPtsR;
nPtsZ = params.nPtsZ;
%might pass these in as params
r = params.r;
z = params.z;
rmat = params.rmat;
%format in vector into nPtsR-2,nPtsZ-2,2 matrix
in = reshape(in,[nPtsR-2,nPtsZ-2,2]);
eta = zeros(nPtsR,nPtsZ);
v = zeros(nPtsR,nPtsZ);
%eta and v at time t
eta(2:end-1,2:end-1) = in(:,:,1);
v(2:end-1,2:end-1) = in(:,:,2);

    %get psi at time t+delta t
    psi = PsiFromEta(eta,params);
    %use this for BCs of v, eta at t
    v = VBCs(v,params);
    
        %step 5

        eta = EtaBCs(eta,psi,params);
     
    %steps 3,4 now get detadt dvdt for t+delta t
    dpsidz = ddz(psi,z);
    dpsidr = ddr(psi,r);
    J = @(x) (dpsidz .* (ddr(x,r))) - (dpsidr .* (ddz(x,z)));
    dvdt = J(v)./rmat + v.*dpsidz./rmat.^2 ;

    temp = eta./rmat;
    %handle the 0/0 problem
    temp(eta==0) =0;
    
    temp(1,:) = temp(2,:);
    detadt =  2*v.*ddz(v,z)./rmat;
    detadt(1,:)= 0;
    detadt = J(temp) +detadt;
    
    out(:,:,1) = detadt(2:end-1,2:end-1);
    out(:,:,2) = dvdt(2:end-1,2:end-1);
    out = out(:);

end


function eta = EtaBCs(eta,psi,params)
r = params.r;
R = params.R;
z = params.z;
etaInit = params.EtaInit;

dr = r(2)-r(1);
    %r=0
    eta(1,:) = etaInit(0,z);
    %r=R
    eta(end,:) = (1/R^2).*((psi(end,:) - psi(end-1,:))/dr) - (1/R) .*((psi(end,:) - 2*psi(end-1,:) + psi(end-2,:))/(dr^2));

    %z=0
    %not enforcing anything
    i = 2:length(r)-1;
    eta(i,1) = ((1./r(i).^2)'.*((psi(i+1,1) - psi(i-1,1))/(2*dr))) ...
        - ((1./r(i))' .*((psi(i+1,1) - 2*psi(i,1) + psi(i-1,1))/(dr^2)));
    %z=Z
    
    eta(:,end) = eta(:,end-1);
end

function v = VBCs(v,params)
r = params.r;
z = params.z;
R = params.R;
VInit = params.VInit;
    %r=0
    %Limit of VInit(r-> 0,z) = 0
    v(1,:) = 0; %VInit(0,z);
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