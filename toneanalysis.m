load tone3.mat
x=audioarray';
x=x(2401:24000);
Sec=2.7;
tau=[0:Fs*Sec-1]/Fs;
dt=0.01;
t=[0:dt:max(tau)];
f=[0:1:500];
addpath C:\Users\MD531\Documents\MATLAB\TFW
sgm=500;
y=abs(Gabor51(x,tau,t,f,sgm));
% sgm=20;
% y2=Gabor51(x,tau,t,f,sgm);
% y=abs(y1).^0.5.*abs(y2).^0.5;
figure(3)
image(t,f,y/max(max(y))*600)   
colormap(gray(256))         
set(gca,'Ydir','normal')    
set(gca,'Fontsize',14)         
xlabel('Time (Sec)','Fontsize',14)             
ylabel('Frequency (Hz)','Fontsize',14)      




