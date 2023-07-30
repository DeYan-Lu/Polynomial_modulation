function [X, t, f,Q1] = Gabor_ub1(x, N, dtau, s, sigma)

% x: input signal
% N: sampling points on f-axis
% dtau: sampling inteval of x
% s: unbalanced sampling parameter
% Q: window size
% sigma: Gaussian standard deviation

T=length(x);
t=[1:s:T]*dtau;                 % output t-axis  %% modified t=[1:s:T-s+1]*dtau;
C=length(t);
df=1/N/dtau;
m=[-fix(N/2):-fix(N/2)+N-1];  %% modified 
f=m*df;   % plot f-axis
B=1.9143/sigma^.5;
Q=min(round(B/dtau),fix((N-1)/2));
xx=[zeros(1,Q),sigma^(0.25)*dtau*x,zeros(1,Q)];
n1=mod(m,N)+1;
X=zeros(N,C);
w=exp(-sigma*pi*([-Q:Q]*dtau).^2);
pha=j*2*pi*m'/N;
Q1=2*Q+1;
for n=1:C
    x1=zeros(1,N);
    x1(1:Q1)=xx((n-1)*s+1:(n-1)*s+Q1).*w;   %% modified  x(n*s:n*s+2*Q)
    X1=fft(x1,N);
    X(:,n)=exp(pha*(Q-(n-1)*s-1)).*X1(n1).';    %% modified  
end
end