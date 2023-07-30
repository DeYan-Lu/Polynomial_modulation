clear; 
close all; clc;
warning('off','all');

%% Functions
addpath('Functions');
addpath('WAV');

%% Input
fprintf('\nAnalyzing input signal......\n');
datanames={'Birds.wav','Cow.wav','Dog.wav','Elephant.wav','Horse.wav','Monkey.wav','Sheep.wav'};
index=6;
filename = datanames{index};
[y,Fs] = wavread(filename);   % true example
x=y(:,1).'; %雙聲道,取第1聲道
% sound(x,Fs);

T=length(x);
N=10000;
dtau=1/Fs;              % sampling interval for input
s_un=fix(Fs*0.01);      % unbalanced sampling parameter
tau=[1:T]*dtau;

% x=zeros(size(tau));             % simple example 
% x(fix(T/7):fix(T/3))=cos(2000*pi*(tau(fix(T/7):fix(T/3))+0.5).^2);
% x((fix(T/2):fix(6*T/7)))=cos(1000*pi*(tau(fix(T/2):fix(6*T/7))-0.8).^3+2000*pi*(tau(fix(T/2):fix(6*T/7))-0.8));

x_ft=fft(x);
x_ft(round(T/2)+1:end)=0;               %把後半部signal(負頻成分切掉)切掉
x_an=ifft(x_ft);                        % analytic signal
x_g=x_an;
j=sqrt(-1);
shift=4;
x_w=x_an.*exp(-j*2*pi*(Fs/shift)*tau);      % translation Fs/4 (頻率軸移動)

tic;
%% Gabor Transform
sigma=2000;
% input:x_w ; N:sampling points on f-axis ; dtau: sampling inteval of x_w
% s_un: unbalanced sampling parameter ; sigma: Gaussian standard deviation
% f1:f-axis ; t:t-axis ;X1:Gabor
[X1, t, f1] = Gabor_ub(x_w, N, dtau, s_un, sigma);

%% Wigner distribution function
% x_w: input signal
% fix(N/2): sampling points on f-axis
% dtau: sampling inteval of x
% s_un: unbalanced sampling parameter
[X2, t, f2] = Wigner_ub(x_w, fix(N/2), dtau, s_un);
%% Axis alignment
x1a=find(abs(f1-f2(1))<Fs/N/2); %x1a位於f1為-5000的位置(第2501行) ;f2(1)=-5000 ;間隔為2
%length(f2)=length(f1)=5000
display(f1)
f1=f1(x1a:x1a+length(f2)-1);
X1=X1(x1a:x1a+length(f2)-1,:);
%要使得f1和f2的軸都是一樣

display(f2(1))
display(x1a)
display(Fs)

figure(1)
C1=1000;
subplot(231)
image(t,f1,abs(X1)/max(max(abs(X1)))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_{xw}(t,f)')

C2=3000;
subplot(232)
image(t,f2,abs(X2)/max(max(abs(X2)))*C2)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Wigner distribution W_{xw}(t,f)')

%% GWT
f3=f2;
alpha=2;
beta=0;
ab=alpha+2*beta;

C3=1e+4;
Cx=abs(X1).^alpha.*abs(X2).^(beta/2);
subplot(233)
image(t,f3,abs(Cx)/max(max(abs(Cx)))*C3)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('GWT C_x(t,f) = |G_{xw}(t,f)|^\alpha |W_{xw}(t,f)|^\beta')

thr_seg=mean(mean(Cx.^(1/ab)*3))^ab;
% constant_K=2e+11;
% thr_seg=(sqrt(sum(x.^2)))^(2*alpha+2*beta)/constant_K;
R_seg=(Cx>=thr_seg);  % 判別式,True=1 ; False=0
% R=Cx1.*(Cx1>=thr_seg);
subplot(234)
image(t,f3,R_seg*255);
colormap(gray(256))
set(gca,'Ydir','normal')
title(['Thresholding R_{seg}(t,f) ( thr = ',num2str(thr_seg),' )'])

set(gcf,'position',[50 100 1000 600]);
toc;

%% Segmentation
fprintf('\nDividing data into parts......\n');
% constant_C=0.8;
% thr_ex=constant_C*N*2/s_un;     % C/df/dt = C/s/dtau/df = C*N*2/s
uc=1/(t(2)-t(1))/(f3(2)-f3(1));     % uncertainty principle
thr_ex=ceil(uc*5); 

tic;

% Dilation operation
r1=1;
r2=6;
[yk,xk]=meshgrid(-r2:r2,-r1:r1);
kernel=((xk/(r1+0.3)).^2+(yk/(r2+0.3)).^2)<1;
[fx,fy]=find(kernel); %找kernel為1的(r,c)座標
center=[r1+1,r2+1];  %求中心點座標
fx=fx-center(1);   %將中心點設為原點(0,0)
fy=fy-center(2);
R_mor=dilation(R_seg,[fx,fy]);

[label,S0]=bwlabel(R_mor); % S0:connected object number = 62
label=label.*R_seg;
display(size(label))
% Sorted by number of labels
count=zeros(1,S0);
for s=1:S0
    count(s)=sum(sum(label==s));
end
display(count)
map_sort=zeros(1,S0+1); %第136行 S0=62 ; S_ex=9
S_ex=sum(count>=thr_ex); %第106行
display(S_ex) %發現有S_ex元素大於等於thr_ex,S_ex=9
[num_label,pos_label]=sort(count,'descend'); % pos_labe描述新排列矩陣num_label在原始矩陣count的位置
display(num_label)
display(pos_label)
map_sort(pos_label(1:S_ex)+1)=[1:S_ex]; %第130行
display(map_sort)
display(size(label))
label_ex=map_sort(label+1);
display(size(label_ex))

subplot(235)
image(t,f3,R_mor*255);
colormap(gray(256))
set(gca,'Ydir','normal')
title('Morphology R_{mor}(t,f)')

subplot(236)
image(t,f3,255*(label_ex~=0))
colormap(gray(256))
set(gca,'Ydir','normal')
title(['Labeled signal ( thr = ',num2str(thr_ex),' )'])

s11=ceil(S_ex^0.5);  s12=ceil(S_ex/s11); %(設定subplot大小)
%了解每個label後的位置
figure(2)
for s=1:S_ex
    subplot(s11,s12,s)   %s11 ,s12分別在第154行
    image(t,f3,255*(label_ex==s))
    colormap(gray(256))
    set(gca,'Ydir','normal')
end
supertitle('Each part of the signal')

toc;

%% Approximation
fprintf('\nApproximating each part of signal......\n');
tic;

% Encoded: S_ex, dtau, Td (t_min_index, t_max_index), Iv (ds_ratio), Coef (x_dk), Pv (polynomial value)
N_order=6;
Td=[]; Iv=[]; Coef=[]; Pv=[];
BW=[]; len_k=[]; len_dk=[]; nmse_div=[];

dt1=0.05;
df1=3/N/dtau;
x_w1=x_w;
nmse_div_thr=0.003;
y_rec=zeros(size(x));
for k=1:S_ex
    R_div=Cx.*(label_ex==k);
    f_ins=freq_ins(R_div,f3); %分母有可能為零
    f_ins_div=f_ins(f_ins~=0);

    % approximation f_ins = a0 + a1*t + a2*t^2 + ..... + aN*t^N
    t_div=t(f_ins~=0);
    a_n=fliplr(polyfit(t_div,f_ins_div,N_order));
    f_poly=polyval(fliplr(a_n),t_div);

    fp=zeros(1,length(t));
    fp(f_ins~=0)=f_poly;
    [fx1,fy1]=find(R_div); %find R_div不為零的地方(依column方向數)
 
    B=max(abs(fp(fy1)-f3(fx1)))+df1;
    BW=[BW,B];
    
%     t_max=max(t_div);
%     t_min=min(t_div);
%     t_max_index=find(tau==t_max);
%     t_min_index=find(tau==t_min);
    
    t_max=min(max(t_div)+dt1,tau(end));
    t_min=max(min(t_div)-dt1,tau(1));
    [~,t_max_index]=min(abs(tau-t_max));
    [~,t_min_index]=min(abs(tau-t_min));
    
    Td=[Td,t_min_index,t_max_index];
    tpoint=tau(t_min_index+round((t_max_index-t_min_index)*[0:N_order]/N_order));
    f_poly1=polyval(a_n,tpoint);  
    Pv=[Pv,f_poly1];

    % Generalized modulation
    t_interval_index=[t_min_index:t_max_index];
    t_interval=tau(t_interval_index);
    fm=[0,a_n./[1:N_order+1]];        % modulation frequency coefficients(經過積分後)
    m_t=exp(j*2*pi*polyval(fliplr(fm),t_interval));
    x_div=x_w1(t_interval_index);
    x_mod=x_div./m_t;
    
    % Cutoff frequency
    X_k_mod=fft(x_mod);
    f_index=round(B/Fs*length(t_interval));
    X_k_mod(f_index+2:end-f_index)=0;
    x_k_mod=ifft(X_k_mod);
    x_k=x_k_mod.*m_t;
    x_w1(t_interval_index)=x_w1(t_interval_index)-x_k;

    min_points=10;
    nmse_d=1;
    while nmse_d>nmse_div_thr  %第179行 = 0.003
        ds_ratio=min(fix(1/2/B/dtau),fix(length(t_interval)/min_points));
        dt=ds_ratio*dtau;
        x_dk=x_k_mod(1:ds_ratio:end);
        
        % reconstruction
        y_k=zeros(size(t_interval));
        n=0;
        for ts=min(t_interval):dt:max(t_interval)
            n=n+1;
            y_dk=x_dk(n)*sinc((t_interval-ts)/dt);
            y_k=y_k+y_dk;
        end
        y_k=y_k.*m_t;
        nmse_d=sum((abs(y_k-x_k)).^2)/sum(abs(x_k).^2);
        min_points=min_points+10;
    end
    y_rec(t_interval_index)=y_rec(t_interval_index)+y_k*2;
    
    Coef=[Coef,x_dk];
    Iv=[Iv,ds_ratio];
    nmse_div=[nmse_div,nmse_d];
    len_k=[len_k,length(x_k)];
    len_dk=[len_dk,length(x_dk)];
    
%     figure()
%     subplot(211)
%     plot(t_interval,real(x_k))
%     xlabel('Time (Sec)')
%     ylabel('Amplitude')
%     title('Signal x_k(t)')
% 
%     subplot(212)
%     plot(t_interval,real(y_k))
%     xlabel('Time (Sec)')
%     ylabel('Amplitude')
%     title('Signal y_k(t)')

end
y_rec=real(y_rec.*exp(j*2*pi*(Fs/4)*tau));
y_rec=y_rec-mean(y_rec)/2;

toc;


%% Result
nmse_total=sum((y_rec-x).^2)/sum(x.^2);

fprintf('\n----------------------------------------------- ');
fprintf('\n< Result >');
fprintf(['\nParts of signal: ',num2str(S_ex)]);
fprintf(['\nTotal NMSE: ',num2str(nmse_total),' = ',num2str(nmse_total*100),' %%']);
fprintf('\n----------------------------------------------- ');
fprintf('\n');

fprintf('\n< Encoded coefficients >');
fprintf('\nS_ex (the number of parts of the signal)');
fprintf('\ndtau (minimum time interval)');
fprintf('\nTd (time interval index for each part, [t_min_index, t_max_index])');
fprintf('\nIv (downsampling ratio, ds_ratio)');
fprintf('\nCoef (downsampled signal for each part, x_dk)');
fprintf('\nPv (polynomial value, encoded a_n)');
fprintf(['\nTotal memory: ',num2str(length([S_ex,dtau,Td,Iv, Coef,Pv]))]);
fprintf('\n');

fprintf('\n< Some infos >');
fprintf('\nNMSE of reconstuction for each part: nmse_div');
fprintf('\nOne-sided bandwidth of each part: BW');
fprintf('\nLength of each part signal: len_k');
fprintf('\nLength of each part downsampling signal: len_dk');
fprintf('\n');

figure(3)
subplot(311)
plot(tau,x);
xlim([0,T*dtau]);
xlabel('Time (Sec)')
ylabel('Amplitude')
title('Original signal x(t)')

subplot(312)
plot(tau,y_rec);
xlim([0,T*dtau]);
xlabel('Time (Sec)')
ylabel('Amplitude')
title('Recovered signal y_{rec}(t)')

subplot(313)
plot(tau,(y_rec-x).^2);
xlim([0,T*dtau]);
xlabel('Time (Sec)')
ylabel('Amplitude')
title(['Error ( y_{rec}(t) - x(t) )^2 ( NMSE = ',num2str(nmse_total),' )'])


[Y, ty, fy] = Gabor_ub(y_rec, N, dtau, s_un, sigma);
Y=abs(Y(fix(N/2)+1:end,:));
fy=fy(fix(N/2)+1:end);

figure(4)
subplot(211)
image(t,f1+Fs/4,abs(X1)/max(max(abs(X1)))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_x(t,f)')

subplot(212)
image(ty,fy,abs(Y)/max(max(abs(Y)))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_{yrec}(t,f)')

%% Sound
% sound(x,Fs);
% sound(y_rec,Fs);

% pause;
