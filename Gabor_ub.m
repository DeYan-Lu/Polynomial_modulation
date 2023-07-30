function [X, t, f] = Gabor_ub(x, N, dtau, s, sigma)

% x: input signal
% N: sampling points on f-axis
% dtau: sampling inteval of x
% s: unbalanced sampling parameter
% Q: window size
% sigma: Gaussian standard deviation

T=length(x);
t=[1:s:T]*dtau;                     % plot t-axis
C=length(t);
df=1/N/dtau;
f=[-fix(N/2):N-fix(N/2)-1]*df;      % plot f-axis
B=1.9143/sqrt(sigma);
Q=min(round(B/dtau),fix((N-1)/2));

xx=[zeros(1,Q),sigma^(0.25)*dtau*x,zeros(1,Q)];
m1=[N-fix(N/2)+1:N,1:N-fix(N/2)];
w=exp(-sigma*pi*([-Q:Q]*dtau).^2);
pha=j*2*pi*[0:N-1]/N;

X=zeros(N,C);
Q1=2*Q+1;
for n=1:C
    x1=zeros(1,N);
    x1(1:Q1)=xx((n-1)*s+1:(n-1)*s+Q1).*w(Q1:-1:1);
    X1=fft(x1,N).*exp(pha*(Q-(n-1)*s+1));
    X(:,n)=X1(m1);
end
end