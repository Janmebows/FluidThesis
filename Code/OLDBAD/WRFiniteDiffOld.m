nPts = 10;
R = 1;
Z = 1;
params.nPts = nPts;
params.R = R;
params.Z = Z;
%obtain the LU factorisation to solve psi with
[L,U,A] = solveAAndDecompose(params);

repR = repmat(linspace(0,R,nPts),1,nPts)';
r = linspace(0,R,nPts);
z = linspace(0,Z,nPts);
[rmat, zmat] = ndgrid(r,z);
%rhs = -rmat.* (sin(2*pi*zmat) .*(-4*pi^2*rmat.^2.*(2*pi-rmat) - 3*rmat));
rhs = -(2*rmat.^2.*(rmat-1) + 2*zmat.*(3*rmat-1).*(zmat-1) - zmat.*(3*rmat-2).*(zmat-1));
%solve L*U*psi = -r*etaV

%r=0
rhs(1,:) = 0;
%r=R
rhs(end,:) = 0;
%z=0
rhs(:,1) =0;
%z=Z
rhs(:,end) =0;
rhsvec = rhs(:);


%A psi = - r eta
% L*y = -r*etaV
y = L\rhsvec;
psiV= U\y;
psi = A\rhsvec;
psi = Vec2Mat(psiV,nPts,nPts);

close all
surf(rmat,zmat,psi)
figure
%surf(rmat,zmat,1 + rmat.^2.*(2*pi-rmat).* sin(2*pi*zmat))
surf(rmat,zmat,-rmat.^2.*(rmat-1).*zmat.*(zmat-1))
xlabel("true sol")
function [L,U,A] = solveAAndDecompose(params)
R = 1;
Z = 1;
nPts = params.nPts;
r = linspace(0,R,nPts);
z = linspace(0,Z,nPts);
dr = r(2) - r(1);
dz = z(2) - z(1);
invdr = 1/dr;
invdz = 1/dz;
A = spalloc(nPts*nPts,nPts*nPts,nPts^2 + 5*(nPts-1)^2);

%loop over j and i
%i can't seem to vectorise this nicely so double loop for now
for j=2:nPts-1
   for i=2:nPts-1
       %b = a + 1
       A(i+(j-1)*nPts,i+1+(j-1)*nPts) = invdr^2 - (invdr./(2*r(i)));
       %b = a
       A(i+(j-1)*nPts,i+(j-1)*nPts) = -2*(invdz^2 + invdr^2);
       %b = a - 1
       A(i+(j-1)*nPts,i-1+(j-1)*nPts) = invdr^2 + (invdr./(2*r(i)));
       %b = a - m
       A(i+(j-1)*nPts,i+(j-2)*nPts) = invdz^2;
       %b = a +m
       A(i+(j-1)*nPts,i+j*nPts) = invdz^2; 
   end    
end
%BCs
for i=1:nPts
    %z = 0 => j=1
    ind = cvtIndex(i,1,nPts);
    A(ind,ind) = 1;
    %z = Z => j = nPts
    ind = cvtIndex(i,nPts,nPts);
    A(ind,ind) = 1;
    %r = 0 => i=1
    ind = cvtIndex(1,i,nPts);
    A(ind,ind) = 1;
    %A(i+((nPts-1)*nPts),i+((nPts-1)*nPts)) = 1;
    %r=R
    ind = cvtIndex(nPts,i,nPts);
    A(ind,ind) = 1;
    %A(ind,ind) = 0.5*R^2;
    
    %A(nPts*i,nPts*i) = 0.5*R^2;
end   
spy(A)
[L, U] = lu(A);
end
% plot3(r,z,eta)
function ind = cvtIndex(i,j,nPts)
    ind = i+(j-1)*nPts;
end

function vec = Mat2Vec(mat)
    vec = mat(:);
end

function mat = Vec2Mat(vec,ndim,mdim)
    mat = zeros(ndim,mdim);
    for j=1:mdim
        lind = cvtIndex(1,j,ndim);
        rind = cvtIndex(0,j+1,mdim);
%         mat(:,j) = vec((j-1)*ndim+1:j*mdim);
        mat(:,j) = vec(lind:rind);
    end
end
