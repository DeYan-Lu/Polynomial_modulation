sigma=100;
tau=2.5;
fs=100;
dt=(2*tau+1)./fs;
x=-tau:dt:tau;
t=-tau:dt:tau;
[eig_function,eig_value]=dpss_func(sigma,tau,x,t);
display(size(eig_function));
display(size(eig_value));

x=1:1:length(eig_value);
y=eig_function(hh,10);
xx=1:0.025:length(eig_value);
yy=spline(x,y,xx);
plot(x,y,'o',xx,yy)

display(eig_value)
u=eig_function(:,3).*eig_function(:,56);
uu=sum(u);
display(uu)

