function [X, t, f] = Wigner_ub(x, N, dtau, s)

% x: input signal
% N: sampling points on f-axis
% dtau: sampling inteval of x
% s: unbalanced sampling parameter
% Q: window size

T=length(x);
t=[1:s:T]*dtau;                     % plot t-axis
C=length(t);
df=1/2/N/dtau;
f=[-fix(N/2):N-fix(N/2)-1]*df;      % plot f-axis

xx=2*dtau*x;
m1=[fix(N/2)+2:N,1:fix(N/2)+1];
pha=j*2*pi*[0:N-1]/N;

X=zeros(N,C);
for n=1:C
    x1=zeros(1,N);
    n1=(n-1)*s+1;
    Q=min(min(n1-1,T-n1),fix((N-1)/2));
    x1(1:2*Q+1)=xx(n1-Q:n1+Q).*conj(xx(n1+Q:-1:n1-Q));
    X1=fft(x1,N).*exp(pha*Q);
    X(:,n)=X1(m1);
end
end