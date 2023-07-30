clear; close all; clc;
%% Input
Fs=250;
t_sec=10;
T=t_sec*Fs;

N=10000;
dt=1/Fs;
B=10;
Q=fix(B/dt);
t=[1:T]*dt;                     % plot t-axis
n1=[fix(N/2)+2:N,1:fix(N/2)+1];

x=exp(j*(t-5).^4-j*5*pi*(t-5).^2);

tic;

%% Gabor Transform
df1=1/N/dt;
f1=[fix(N/2)+1-N:fix(N/2)]*df1;   % plot f-axis
sigma=15;
xx=[zeros(1,Q),sigma^(0.25)*dt*x,zeros(1,Q)];

X1=zeros(T,N);
m=[1:N];
w=exp(-sigma*pi*([-Q:Q]*dt).^2);

for n=1:T
    x1=zeros(1,N);
    x1(1:2*Q+1)=xx(n:n+2*Q).*w(2*Q+1:-1:1);
    X=fft(x1);
    X(m)=exp(i*2*pi*(Q-n+1)*(m-1)/N).*X(m);
    X1(n,:)=X(n1);
end

figure()
C1=400;
subplot(221)
image(t,f1,abs(X1.')/max(max(abs(X1.')))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_x(t,f)')

%% Wigner distribution function
df2=1/2/N/dt;
f2=[fix(N/2)+1-N:fix(N/2)]*df2;
Q=max(fix(T/2),fix((N-1)/2));
xx=[zeros(1,Q),2*dt*x,zeros(1,Q)];

X2=zeros(T,N);
m=[1:N];

for n=1:T
    x2=zeros(1,N);
    Q1=min(min(n-1,T-n),fix((N-1)/2));
    x2(1:2*Q1+1)=xx(n+Q-Q1:n+Q+Q1).*conj(xx(n+Q+Q1:-1:n+Q-Q1));
    X=fft(x2,N);
    X(m)=exp(i*2*pi*Q1*(m-1)/N).*X(m);
    X2(n,:)=X(n1);
end

C2=2000;
subplot(222)
image(t,f2,abs(X2.')/max(max(abs(X2.')))*C2)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Wigner distribution W_x(t,f)')

%% GWT
X1t=zeros(size(X1));
X1t(:,1:2:end)=X1(:,fix(N/4)+1:round(3*N/4));
X1t(:,2:2:end)=X1(:,fix(N/4)+1:round(3*N/4));
sp=abs(X1t).^2;

C3=2500;
Cx1=X1t.*X2;
subplot(223)
image(t,f2,abs(Cx1.')/max(max(abs(Cx1.')))*C3)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('GWT C_{x1}(t,f) = G_x(t,f) W_x(t,f)')

C4=2500;
alpha=2.5;
beta=0.7;
thr=1e-4;
Cx2=(sp.^alpha).*(sp.^alpha>thr).*abs(X2.^beta);
subplot(224)
image(t,f2,abs(Cx2.')/max(max(abs(Cx2.')))*C4)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('GWT C_{x2}(t,f) = thr[SP_x^\alpha(t,f)] |W_x^\beta(t,f)|')

set(gcf,'position',[150 100 700 600]);
toc;

%% Generalized modulation
f_ins=freq_ins(X1.', f1);

% approximation f_ins = a0 + a1*t + a2*t^2 + ..... + aN*t^N
N_order=5;

% approximation by matrix computation
A=zeros(N_order+1,N_order+1);
b=zeros(N_order+1,1);
for i=1:N_order+1
    A(i,:)=sum(repmat(t.',1,N_order+1).^repmat((i-1):(i+N_order-1),T,1));
    b(i,:)=sum(f_ins.*(t.^(i-1)));
end
a_n=A\b;
y=polyval(fliplr(a_n(:).'),t);

% % approximation by polyfit
% a_n=fliplr(polyfit(t,f_ins,N_order)).'
% y=polyval(fliplr(a_n(:).'),t);

figure()
subplot(221)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( G_x(t,f) )')
legend('f_{ins}(t)','approximation curve')

f_ins=freq_ins(X2.', f2);
a_n=fliplr(polyfit(t,f_ins,N_order)).';
y=polyval(fliplr(a_n(:).'),t);

subplot(222)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( W_x(t,f) )')
legend('f_{ins}(t)','approximation curve')

f_ins=freq_ins(Cx1.', f2);
a_n=fliplr(polyfit(t,f_ins,N_order)).';
y=polyval(fliplr(a_n(:).'),t);

subplot(223)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( C_{x1}(t,f) )')
legend('f_{ins}(t)','approximation curve')

f_ins=freq_ins(Cx2.', f2);
a_n=fliplr(polyfit(t,f_ins,N_order)).';
y=polyval(fliplr(a_n(:).'),t);

subplot(224)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( C_{x2}(t,f) )')
legend('f_{ins}(t)','approximation curve')

set(gcf,'position',[250 100 700 600]);
