clear; close all; clc;
%% Image
filename='peppers256.bmp';
% filename1='MaxAlpha=2.bmp';
im_orig=double(imread(filename))/255;
im_new=zeros(size(im_orig));
im1=zeros(size(im_orig));
im2=zeros(size(im_orig));
[M1,N1]=size(im_orig);
im=[im_orig im_orig(:,N1)];

%% Distortion
beta=0.7;
alpha_min=0.2;
m=[1:M1]';
% alpha=alpha_min.^(((M1-m)/(M1-1)).^beta);
alpha=1+alpha_min-alpha_min.^(((m-1)/(M1-1)).^beta);

for m=1:M1
    N2=fix(N1*alpha(m));
    N2_trans=fix((N1-N2)/2);
    for n=1:N2
        nn=fix(n/alpha(m));
        x=n/alpha(m)-nn;
        im_new(m,n+N2_trans)=(1-x)*im(m,nn)+x*im(m,nn+1);
    end
end

%% Recovery
% Linear interpolation
for m=1:M1
    index1=find(im_new(m,:),1);
    index2=find(im_new(m,:),1,'last');
    temp=[im_new(m,index1) im_new(m,index1:index2) im_new(m,index2)];
    length_row=index2-index1+1;
    alpha_re=length_row/N1;
    for n=1:N1
        nn=fix(n*alpha_re)+1;
        x=n*alpha_re-nn+1;
        im1(m,n)=(1-x)*temp(nn)+x*temp(nn+1);
    end
end

% Spline
for m=1:M1
    index1=find(im_new(m,:),1);
    index2=find(im_new(m,:),1,'last');
    temp=im_new(m,index1:index2);
    length_row=index2-index1+1;
    im2(m,:)=spline(linspace(0,1,length_row),temp,linspace(0,1,N1));
end

%% Result
NRMSE1=sqrt(sum(sum((im1-im_orig).^2))/sum(sum(im_orig.^2)));
NRMSE2=sqrt(sum(sum((im2-im_orig).^2))/sum(sum(im_orig.^2)));
SSIM1=SSIM(im1,im_orig,0.00005,0.00005);
SSIM2=SSIM(im2,im_orig,0.00005,0.00005);

%% Figure
figure
imagesc(im_orig);
colormap(gray(256));
axis image
title('Original Image');
figure
imagesc(im_new);
colormap(gray(256));
axis image
title('Distorted Image');

figure
imagesc(im1);
colormap(gray(256));
axis image
title({['Recovery using linear interpolation by rows'];['NRMSE=',num2str(NRMSE1),', SSIM=',num2str(SSIM1)]});
figure
imagesc(im2);
colormap(gray(256));
axis image
title({['Recovery using spline function by rows'];['NRMSE=',num2str(NRMSE2),', SSIM=',num2str(SSIM2)]});
