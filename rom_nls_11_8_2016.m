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
figure(2)
plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2])

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

% sampling locations P 
P = eye(n)

% P = zeros(10,n); % spits out 10 x 3 matrix at the backend 
% for j = 1:10
%    P(j,128+(2*j-1)) = 1;
% end
% 
% figure(31),spy(P)
 PNL = P*phi; % 10 x 3 matrix 
 figure(33),spy(PNL)

[t,asol]=ode45('a_rhs_11_8',t,a0,[],phi,Lr,PNL);
% asol has columns, a1, a2, a3 , in phi1 phi2 phi3 direction 
  
us=zeros(n,length(t));
for j=1:length(t)
    for jj=1:r
       us(:,j)=us(:,j) + asol(j,jj)*phi(:,jj);  % r-rank reconstruction
    end
end
figure(4)
surfl(x,t,abs(us.')); shading interp, colormap(hot);








