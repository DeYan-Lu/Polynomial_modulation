function [eig_function,eig_value] = dpss_func(sigma,tau,x,t)
dt=x(2)-x(1);
Q=round(tau/dt);
p=x./dt;
q=t./dt;
%parameter:
R=sigma/pi;
A=zeros(length(q)+1,length(p)+1);
for xx=1:1:2*Q+1
    for yy=1:1:2*Q+1
        A(xx,yy)=A(xx,yy)+R.*sinc(R*(xx-yy)*dt);
    end
end
[eig_function,eig_value_1]=eig(A);
eig_value_2=zeros(1,length(q));

for zz=2*Q+1:-1:1
    eig_value_2(1,2*Q+1-zz+1)=eig_value_1(zz,zz);
end
% normalized eigenvalue
eig_max=max(eig_value_2);
eig_value=eig_value_2./eig_max;
end
