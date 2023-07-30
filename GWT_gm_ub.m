clear; close all; clc;

%% Functions
addpath('Functions');

%% Input
Fs=250;
t_sec=10;
T=t_sec*Fs;

N=10000;
dtau=1/Fs;      % sampling interval for input
s=5;            % unbalanced sampling parameter
tau=[1:T]*dtau;

x=exp(j*(tau-5).^4-j*5*pi*(tau-5).^2);

tic;

%% Gabor Transform
sigma=20;
B=1;            % window size
Q=fix(B/dtau);

[X1, t, f1] = Gabor_ub(x, N, dtau, s, Q, sigma);

figure()
C1=400;
subplot(221)
image(t,f1,abs(X1)/max(max(abs(X1)))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_x(t,f)')

%% Wigner distribution function

[X2, t, f2] = Wigner_ub(x, N, dtau, s);

C2=2000;
subplot(222)
image(t,f2,abs(X2)/max(max(abs(X2)))*C2)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Wigner distribution W_x(t,f)')

%% GWT
X1t=zeros(size(X1));
X1t(1:2:end,:)=X1(fix(N/4)+1:round(3*N/4),:);
X1t(2:2:end,:)=X1(fix(N/4)+1:round(3*N/4),:);
sp=abs(X1t).^2;

C3=2500;
Cx1=X1t.*X2;
subplot(223)
image(t,f2,abs(Cx1)/max(max(abs(Cx1)))*C3)
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
image(t,f2,abs(Cx2)/max(max(abs(Cx2)))*C4)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('GWT C_{x2}(t,f) = thr[SP_x^\alpha(t,f)] |W_x^\beta(t,f)|')

set(gcf,'position',[150 100 700 600]);
toc;

%% Generalized modulation
f_ins=freq_ins(X1,f1);

% approximation f_ins = a0 + a1*t + a2*t^2 + ..... + aN*t^N
N_order=5;

% approximation by matrix computation
A=zeros(N_order+1,N_order+1);
b=zeros(N_order+1,1);
for i=1:N_order+1
    A(i,:)=sum(repmat(t.',1,N_order+1).^repmat((i-1):(i+N_order-1),length(t),1));
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

f_ins=freq_ins(X2, f2);
a_n=fliplr(polyfit(t,f_ins,N_order)).';
y=polyval(fliplr(a_n(:).'),t);

subplot(222)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( W_x(t,f) )')
legend('f_{ins}(t)','approximation curve')

f_ins=freq_ins(Cx1, f2);
a_n=fliplr(polyfit(t,f_ins,N_order)).';
y=polyval(fliplr(a_n(:).'),t);

subplot(223)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( C_{x1}(t,f) )')
legend('f_{ins}(t)','approximation curve')

f_ins=freq_ins(Cx2, f2);
a_n=fliplr(polyfit(t,f_ins,N_order)).';
y=polyval(fliplr(a_n(:).'),t);

subplot(224)
plot(t,f_ins,t,y)
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Approximation of f_{ins}(t) ( C_{x2}(t,f) )')
legend('f_{ins}(t)','approximation curve')

set(gcf,'position',[250 100 700 600]);

%% Modulation
fm=[0,a_n.'./[1:N_order+1]];        % modulation frequency coefficients
m_t=exp(-j*2*pi*polyval(fliplr(fm),tau));
y=x.*m_t;

[X, t, f] = Gabor_ub(y, N, dtau, s, Q, sigma);
figure()
image(t,f,abs(X)/max(max(abs(X)))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_y(t,f) ( y(t) = m(t)x(t) )')
