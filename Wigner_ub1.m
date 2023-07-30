function [X, t, f] = Wigner_ub1(x, N, dtau, s)

% x: input signal
% N: sampling points on f-axis
% dtau: sampling inteval of x
% s: unbalanced sampling parameter
% Q: window size

T=length(x);
t=[1:s:T]*dtau;                 % plot t-axis
C=length(t);
df=1/2/N/dtau;
m=[-fix(N/2):-fix(N/2)+N-1];  %% modified 
f=m*df;   % plot f-axis
xx=[2*dtau*x];  %% modified
m1=mod(m,N)+1;
X=zeros(N,C);
pha=j*2*pi*m'/N;
ns=round(t/dtau);
N0=fix((N-1)/2);
for n=1:C
    x1=zeros(1,N);
    n1=ns(n);
    Q1=min(min(n1-1,T-n1),N0); %% modified 
    x1(1:2*Q1+1)=xx(n1-Q1:n1+Q1).*conj(xx(n1+Q1:-1:n1-Q1)); %% modified 
    X1=fft(x1,N);
    X(:,n)=exp(pha*Q1).*X1(m1).';
end
end