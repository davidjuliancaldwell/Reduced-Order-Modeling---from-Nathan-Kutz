%% lecture 1
clear all;close all;clc

L=20; n = 80;
x2 = linspace(-L/2,L/2,n+1); 
x=x2(1:n);

sigma = 1;
u=exp(-x.^2);
ut = fft(u); 
k = (2*pi/L)*[0:(n/2-1) -n/2:-1];

subplot(2,1,1)
plot(x,u,'ro-','Linewidth',[2])
subplot(2,1,2)
plot(fftshift(k),fftshift(abs(ut)),'bo-')

%% lecture 2

clear all;close all;clc

L=30; n=512; 
x2 = linspace(-L/2,L/2,n+1);
x=x2(1:n);
t=0:0.2:20;
V=(x.^2).'; % x^2 potential, as a column 

k=(2*pi/L)*[0:n/2-1 -n/2:-1].';

u=exp(-0.2*x.^2);
ut=fft(u);

[t,utsol] = ode45('harm_rhs',t,ut,[],k,V);

for j = 1:length(t)
    usol(j,:)=ifft(utsol(j,:));
end
%%
surf(x,t,abs(usol)),shading interp

X = usol.';
[u,s,v] = svd(X,'econ');

% plot singular values
figure(2),plot(diag(s),'ko','Linewidth',[2])

% only even dynamics! POD is limited by sampling strategy 
figure(3),plot(abs(u(:,1:3)))