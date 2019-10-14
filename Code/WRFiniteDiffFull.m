%initialise eta and v
eta = 
v = 




dpsidz = ddz(psi,z);
dpsidr = ddr(psi,r);
J = @(x) dpsidz .* (ddz(x,z)) - dpsidr .* (ddr(x,r));
dvdt = J(v)./r + v.*dpsidz./r.^2 ;
detadt = J(eta./r) + 2*v.*dvdz./r;



nPtsR = 500;
nPtsZ = 500;
R = 1;
Z = 1;
params.nPtsR = nPtsR;
params.nPtsZ = nPtsZ;
params.R = R;
params.Z = Z;
%precompute LU factors

[L,U,A] = solveAAndDecompose(params);

repR = repmat(linspace(0,R,nPtsR),1,nPtsZ)';
r = linspace(0,R,nPtsR);
z = linspace(0,Z,nPtsZ);
[rmat, zmat] = ndgrid(r,z);
rhs = -(2*rmat.^2.*(rmat-1) + 2*zmat.*(3*rmat-1).*(zmat-1) - zmat.*(3*rmat-2).*(zmat-1));


%%BCs----
%r=0
rhs(1,:) = 0;
%r=R
rhs(end,:) = 0;
%z=0
rhs(:,1) =0;
%z=Z
rhs(:,end) =0;

%solve L*U*psi = -r*etaV
rhsvec = rhs(:);

y = L\rhsvec;
psiV= U\y;
psi = Vec2Mat(psiV,nPtsR,nPtsZ);

close all
surf(rmat,zmat,psi)
figure
truSol = -rmat.^2.*(rmat-1).*zmat.*(zmat-1);
surf(rmat,zmat,truSol)
xlabel("true sol")

% diff = abs(psi - truSol);
% figure
% surf(rmat,zmat,diff)
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
spy(A)
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
    %A(i+((nPts-1)*nPts),i+((nPts-1)*nPts)) = 1;
    %r=R
    ind = cvtIndex(nPtsR,j,nPtsR);
    A(ind,ind) = 1;
    %A(ind,ind) = 0.5*R^2;
    
    %A(nPts*i,nPts*i) = 0.5*R^2;
end   
spy(A)
[L, U] = lu(A);
end
% plot3(r,z,eta)
function ind = cvtIndex(i,j,nPtsR)
    ind = i+(j-1)*nPtsR;
end

function vec = Mat2Vec(mat)
    vec = mat(:);
end

function mat = Vec2Mat(vec,ndim,mdim)
    mat = zeros(ndim,mdim);
    for j=1:mdim
        lind = cvtIndex(1,j,ndim);
        rind = cvtIndex(0,j+1,ndim);
%         mat(:,j) = vec((j-1)*ndim+1:j*mdim);
        mat(:,j) = vec(lind:rind);
    end
end




function dxdz = ddz(x,z)
    dz = z(2) - z(1);
    j = 2:length(z)-1;
    dxdz = (x(:,j+1) - x(:,j-1))/(2*dz);
end

function dxdr = ddr(x,r)
    dr = r(2) - r(1);
    i = 2:length(r)-1;
    dxdr = (x(i+1,:) - x(i-1,:))/(2*dr);
end