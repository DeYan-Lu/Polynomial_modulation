clear; 
close all; clc;
warning('off','all');

%% Functions
addpath('Functions');
addpath('WAV');

%% Input
fprintf('\nAnalyzing input signal......\n');
datanames={'Birds.wav','Cow.wav','Dog.wav','Elephant.wav','Horse.wav','Monkey.wav','Sheep.wav'};
index=1;
filename = datanames{index};
[y,Fs] = wavread(filename);   % true example
xa=y(:,1).';            % input signal 
% sound(x,Fs);
T=length(xa);
xf=fft(xa);
xf(floor(xf/2)+2:T)=0;  
x1=ifft(xf);    % analytic signal


N=10000;
dtau=1/Fs;      % sampling interval for input
s_un=50;        % unbalanced sampling parameter
tau=[1:T]*dtau;

% x=zeros(size(tau));             % simple example 
% x(fix(T/7):fix(T/3))=cos(2000*pi*(tau(fix(T/7):fix(T/3))+0.5).^2);
% x((fix(T/2):fix(6*T/7)))=cos(1000*pi*(tau(fix(T/2):fix(6*T/7))-0.8).^3+2000*pi*(tau(fix(T/2):fix(6*T/7))-0.8));

figure(1)
subplot(221)
plot(tau,x);
xlim([0,T*dtau]);
xlabel('Time (Sec)')
ylabel('Amplitude')
title('Signal x(t)')

%% Gabor Transform
tic;

sigma=3000;
B_window=0.1;           % window size
Q=fix(B_window/dtau);

[X1, t, f1] = Gabor_ub(x, N, dtau, s_un, Q, sigma);

C1=1000;
subplot(222)
image(t,f1,abs(X1)/max(max(abs(X1)))*C1)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title('Gabor transform G_x(t,f)')

%% Wigner distribution function

% [X2, t, f2] = Wigner_ub(x, N, dtau, s_un);
% 
% C2=1000;
% subplot(223)
% image(t,f2,abs(X2)/max(max(abs(X2)))*C2)
% colormap(gray(256))
% set(gca,'Ydir','normal')
% xlabel('Time (Sec)')
% ylabel('Frequency (Hz)')
% title('Wigner distribution W_x(t,f)')

% %% GWT
% X1t=zeros(fix(N/2),size(X1,2));
% n1=[round(N/2)+1:round(N/2)+round(N/4)];
% n2=[round(N/2)+1:round(N/2)+fix(N/4)];
% X1t(1:2:end,:)=X1(n1,:);
% X1t(2:2:end,:)=X1(n2,:);
% X2t=X2(round(N/2)+1:end,:);
% f3=f2(round(N/2)+1:end);
% sp=abs(X1t).^2;
% 
% C3=1e+4;
% Cx1=X1t.*X2t;
% subplot(234)
% image(t,f3,abs(Cx1)/max(max(abs(Cx1)))*C3)
% colormap(gray(256))
% set(gca,'Ydir','normal')
% xlabel('Time (Sec)')
% ylabel('Frequency (Hz)')
% title('GWT C_{x1}(t,f) = G_x(t,f) W_x(t,f)')
% 
% C4=1e+8;
% alpha=2.5;
% beta=0.7;
% Cx2=(sp.^alpha).*abs(X2t.^beta);
% subplot(235)
% image(t,f3,abs(Cx2)/max(max(abs(Cx2)))*C4)
% colormap(gray(256))
% set(gca,'Ydir','normal')
% xlabel('Time (Sec)')
% ylabel('Frequency (Hz)')
% title('GWT C_{x2}(t,f) = SP_x^\alpha(t,f) |W_x^\beta(t,f)|')

% f3=f1(round(N/2)+1:end);
% Cx2=abs(X1(round(N/2)+1:end,:).^2);
% C4=C1;
% constant_K=1e+15;
% thr_seg=(sqrt(sum(x.^2)))^(2*alpha)/constant_K;

% constant_K=1e+31;
% constant_K=1e+28;
% thr_seg=(sqrt(sum(x.^2)))^(2*alpha+2*beta)/constant_K;
df0=2;  dfs=round(df0/(f1(2)-f1(1)));
f3a=find(f1>=0);  f3a=f3a(1:dfs:end);
f3=f1(f3a);
X1a=abs(X1(f3a,:));
thr_seg=mean(mean(X1a))*0.95;
R=X1a.*(X1a>=thr_seg);
subplot(223)
image(t,f3,(R~=0)*256)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title(['Thresholding R(t,f) ( thr = ',num2str(thr_seg),' )'])
subplot(224)
image(t,f3,R/max(max(R))*1000)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title(['Thresholding R(t,f) ( thr = ',num2str(thr_seg),' )'])

set(gcf,'position',[50 100 1000 600]);

toc;
return
%% Segmentation
fprintf('\nDividing data into parts......\n');
constant_C=0.8;
[sz1,sz2]=size(X1a);
thr_ex= 5/((f3(2)-f3(1))*s_un/Fs); % from the uncertainty principle
% thr_ex=constant_C*N*2/s_un;     % C/df/dt = C/s/dtau/df = C*N*2/s
tic;
[label,S0]=bwlabel(R);
mask_ex=zeros(size(R));
R_ex=mask_ex;
S1=0;
for s=1:S0
    bw=(label==s);
    count=sum(sum(bw));
    if count>=thr_ex        
        S1=S1+1;
        mask_ex=mask_ex+bw;
        R_ex=R_ex+S1*bw;
    end
end
% R_ex=R.*mask_ex;
% [~,S1]=bwlabel(mask_ex);
figure(2)
subplot(221)
image(mask_ex*255)
colormap(gray(256))

% kernel=[     -1,-2;-1,-1;-1,0;-1,1;-1,2;
%         0,-3; 0,-2; 0,-1; 0,0; 0,1; 0,2; 0,3;
%               1,-2; 1,-1; 1,0; 1,1; 1,2;        ];
kernel=[     -1,-1;-1,0;-1,1;
        0,-2; 0,-1; 0,0; 0,1; 0,2;
              1,-1; 1,0; 1,1;       ];
% kernel=[     -1,0;
%         0,-1; 0,0; 0,1;
%               1,0;       ];
return
t_mor=2;
mask_mor=closing(mask_ex,kernel,t_mor);
R_mor=R.*mask_mor;

[label,S2]=bwlabel(mask_mor);

figure()
subplot(211)
% image(t,f3,mask_ex*256)
image(t,f3,R_ex/max(max(R_ex))*C4)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title(['R_{ex}(t,f) ( excluding small area of R(t,f) , thr = ',num2str(thr_ex),' )'])

subplot(212)
image(t,f3,mask_mor*256)
% image(t,f3,R_mor/max(max(R_mor))*C4)
colormap(gray(256))
set(gca,'Ydir','normal')
xlabel('Time (Sec)')
ylabel('Frequency (Hz)')
title(['mask_{mor}(t,f) ( closing of mask_{ex}(t,f) )'])

toc;

%% Approximation
fprintf('\nApproximating each part of signal......\n');
tic;

nmse_div_thr=0.003;
BW=[];
len_k=[];
len_dk=[];
nmse_div_his=[];
x_trun=zeros(size(x));
x_cut=zeros(size(x));
y_rec=zeros(size(x));
for k=1:S2
    R_div=R.*(label==k);
    f_ins=freq_ins(R_div,f3);
    f_ins_div=f_ins(f_ins~=0);

    % approximation f_ins = a0 + a1*t + a2*t^2 + ..... + aN*t^N
    N_order=6;
    t_div=t(f_ins~=0);
    a_n=fliplr(polyfit(t_div,f_ins_div,N_order)).';
    f_poly=polyval(fliplr(a_n(:).'),t_div);

    ff=repmat(f3(:),1,size(R_div,2)).*(R_div~=0);
    ff(ff==0)=NaN;
    f_max=max(ff);
    f_max=f_max(~isnan(f_max));
    f_min=min(ff);
    f_min=f_min(~isnan(f_min));
    B=max(max(abs(f_max-f_poly),abs(f_min-f_poly)));
    BW=[BW,B];
    
    t_max=max(t_div);
    t_min=min(t_div);
    t_max_index=find(tau==t_max);
    t_min_index=find(tau==t_min);

%     figure()
%     subplot(211)
%     image(t,f3,R_div/max(max(R_div))*C4)
%     colormap(gray(256))
%     set(gca,'Ydir','normal')
%     xlabel('Time (Sec)')
%     ylabel('Frequency (Hz)')
%     title(['R_{div}(t,f) ( Division of R_{mor}(t,f) )'])
% 
%     subplot(212)
%     plot(t_div,f_ins_div,t_div,f_poly,t_div,f_max,t_div,f_min)
%     xlabel('Time (Sec)')
%     ylabel('Frequency (Hz)')
%     legend('f_{ins}(t)','Approx.','f_{max}(t)','f_{min}(t)')

    % Generalized modulation
    t_interval=[t_min:dtau:t_max];
    t_interval_index=[t_min_index:t_max_index];
%     zp=fix(length(t_interval)*0.1);%
%     t_interval=[zeros(1,zp),t_min:dtau:t_max,zeros(1,zp)];%
%     t_interval_index=[t_min_index-zp:t_max_index+zp];%
    fm=[0,a_n.'./[1:N_order+1]];        % modulation frequency coefficients
    m_t=exp(j*2*pi*polyval(fliplr(fm),t_interval));
    x_div=x(t_interval_index);
    x_mod=x_div./m_t;
    x_trun(t_interval_index)=x_div;
    
    X_k_mod=fft(x_mod);
    f_index=fix(B/Fs*length(t_interval))+1;
    X_k_mod(f_index+1:end-f_index)=0;
    x_k_mod=ifft(X_k_mod);
    x_k=x_k_mod.*m_t;
%     x_k=x_k(zp+1:end-zp);%
%     t_interval=[t_min:dtau:t_max];%
%     t_interval_index=[t_min_index:t_max_index];%
%     m_t=exp(j*2*pi*polyval(fliplr(fm),t_interval));%
    x_cut(t_interval_index)=x_cut(t_interval_index)+x_k*2;

    min_points=10;
    nmse_div=1;
    while nmse_div>nmse_div_thr
        ds_ratio=min(fix(1/2/B/dtau),fix(length(t_interval)/min_points));
        dt=ds_ratio*dtau;


        x_dk=x_k_mod(1:ds_ratio:end);
        y_k=zeros(size(t_interval));
        n=0;
        for ts=t_min:dt:t_max
            n=n+1;
            y_dk=x_dk(n)*sinc((t_interval-ts)/dt);
            y_k=y_k+y_dk;
        end
        y_k=y_k.*m_t;
        nmse_div=sum((real(y_k)-real(x_k)).^2)/sum(real(x_k).^2);
        min_points=min_points+10;
    end
    y_rec(t_interval_index)=y_rec(t_interval_index)+y_k*2;
    len_k=[len_k,length(x_k)];
    len_dk=[len_dk,length(x_dk)];
    nmse_div_his=[nmse_div_his,nmse_div];
    
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
x_cut=real(x_cut);
y_rec=real(y_rec);

toc;


%% Result
nmse_trun=sum((x_trun-x).^2)/sum(x.^2);
nmse_cut=sum((x_cut-x_trun).^2)/sum(x_trun.^2);
nmse_rec=sum((y_rec-x_cut).^2)/sum(x_cut.^2);
nmse_total=sum((y_rec-x).^2)/sum(x.^2);
ssim=SSIM(y_rec,x,0.00005,0.00005);

fprintf('\n----------------------------------------------- ');
fprintf('\n< Result >');
fprintf(['\nParts of signal: ',num2str(S2)]);
fprintf(['\nNMSE of truncating the signal: ',num2str(nmse_trun)]);
fprintf(['\nNMSE of cutoff on frequency: ',num2str(nmse_cut)]);
fprintf(['\nNMSE of reconstucting the signal: ',num2str(nmse_rec)]);
fprintf(['\nTotal NMSE: ',num2str(nmse_total)]);
% fprintf(['\nTotal SSIM: ',num2str(ssim)]);
fprintf('\n----------------------------------------------- ');
fprintf('\n');

fprintf('\n< Some infos >');
fprintf('\nNMSE of reconstuction for each part: nmse_div_his');
fprintf('\nOne-sided bandwidth of each part: BW');
fprintf('\nLength of each part signal: len_k');
fprintf('\nLength of each part downsampling signal: len_dk');
fprintf('\n');

figure()
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
% plot(tau,(y_rec-x).^2./(x.^2));
xlim([0,T*dtau]);
xlabel('Time (Sec)')
ylabel('Amplitude')
title(['Error ( y_{rec}(t) - x(t) )^2 ( NMSE = ',num2str(nmse_total),' )'])


[Y, ty, fy] = Gabor_ub(y_rec, N, dtau, s_un, Q, sigma);

figure()
subplot(211)
image(t,f1,abs(X1)/max(max(abs(X1)))*C1)
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
%%
% clc;
% k=1;
% R_div=R.*(label==k);
% f_ins=freq_ins(R_div,f3);
% f_ins_div=f_ins(f_ins~=0);
% 
% N_order=6;
% t_div=t(f_ins~=0);
% a_n=fliplr(polyfit(t_div,f_ins_div,N_order)).';
% f_poly=polyval(fliplr(a_n(:).'),t_div);
% 
% ff=repmat(f3(:),1,size(R_div,2)).*(R_div~=0);
% ff(ff==0)=NaN;
% f_max=max(ff);
% f_max=f_max(~isnan(f_max));
% f_min=min(ff);
% f_min=f_min(~isnan(f_min));
% B=max(max(abs(f_max-f_poly),abs(f_min-f_poly)))
% 
% t_max=max(t_div);
% t_min=min(t_div);
% t_max_index=find(tau==t_max);
% t_min_index=find(tau==t_min);
% 
% t_interval=[t_min:dtau:t_max];
% t_interval_index=[t_min_index:t_max_index];
% fm=[0,a_n.'./[1:N_order+1]];        % modulation frequency coefficients
% m_t=exp(j*2*pi*polyval(fliplr(fm),t_interval));
% x_div=x(t_interval_index);
% x_mod=x_div./m_t;
% 
% X_k_mod=fft(x_mod);
% f_index=fix(B/Fs*length(t_interval))+1
% 
% aa=abs(X_k_mod);
% figure()
% plot(aa)