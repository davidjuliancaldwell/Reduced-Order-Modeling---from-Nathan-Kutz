clear all; close all; clc

L=30; n=256;
x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
t=linspace(0,2*pi,101);

u=2*sech(x); % initial conditions 
ut=fft(u);

% solving n differential equations simultaneously, 256 diff eqs 
[t,utsol]=ode45('nls_rhs',t,ut,[],k);

for j=1:length(t)
   usol(j,:)=ifft(utsol(j,:));
end

surfl(x,t,abs(usol)); shading interp, colormap(hot);


%% SVD 

X=usol.';
[u,s,v]=svd(X);
[u2,s2,v2] = svd((abs(X).^2).*X);
figure(2)
plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2])

figure(22)
plot(diag(s2)/sum(diag(s2)),'ko','Linewidth',[2])
figure(44)
plot(real(u2(:,1:3)),'Linewidth',[2])

figure(3)
plot(real(u(:,1:3)),'Linewidth',[2])

figure(5)
plot(real(v(:,1:3)),'Linewidth',[2])

r=3;  % rank of projection
phi=u(:,1:r);  % Phi_r POD modes

for j=1:r
  phixx(:,j)=  -ifft((k.^2).*fft(phi(:,j))); % second derivatives % this is the L Phi term % stack em up all up 
  a0(j)= 2*sech(x)*conj(phi(:,j));  % projection of initial conditions
end
Lr= (i/2)* phi'*phixx;  % Low-rank approximation of linear term % this is the 2x2 that you only have to do once, and then feed it in 

%DIEM ALGORITHM
phiNL=u2(:,1:4);
% step 1
[mx,nmax] = max(phiNL(:,1)); % get max value and index 
P=zeros(n,1); P(nmax)=1;

% step 2
c = (P.'*phiNL(:,1))\(P.'*phiNL(:,2));
R = phiNL(:,2)-phiNL(:,1)*c;
[mx,nmax]=max(R);
    P2 = zeros(n,1); P2(nmax)=1;
    P = [P P2];

% step 3 - first two columns , backslash to solve 
c = (P.'*phiNL(:,1:2))\(P.'*phiNL(:,3));
R = phiNL(:,3)-phiNL(:,1:2)*c;
[mx,nmax]=max(R);
    P2 = zeros(n,1); P2(nmax)=1;
    P = [P P2];

% step 4
c = (P.'*phiNL(:,1:3))\(P.'*phiNL(:,4));
R = phiNL(:,4)-phiNL(:,1:3)*c;
[mx,nmax]=max(R);
    P2 = zeros(n,1); P2(nmax)=1;
    P = [P P2];

 spy(P)
 
 Piv = phiNL*inv(P.'*phiNL); 
 
 % multiply it by phi 
 
 P

    % sampling locations P 
%P = eye(n);

% P = zeros(10,n); % spits out 10 x 3 matrix at the backend 
% for j = 1:10
%    P(j,128+(2*j-1)) = 1;
% end
% 
% figure(31),spy(P)
 PNL = P*phi; % 10 x 3 matrix 
 figure(33),spy(PNL)

[t,asol]=ode45('a_rhs_11_23',t,a0,[],phi,Lr,PNL);
% asol has columns, a1, a2, a3 , in phi1 phi2 phi3 direction 
  
us=zeros(n,length(t));
for j=1:length(t)
    for jj=1:r
       us(:,j)=us(:,j) + asol(j,jj)*phi(:,jj);  % r-rank reconstruction
    end
end
figure(4)
surfl(x,t,abs(us.')); shading interp, colormap(hot);








