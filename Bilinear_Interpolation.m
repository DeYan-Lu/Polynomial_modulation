clear; close all; clc;
%% Image
filename='peppers256.bmp';
im_orig=double(imread(filename))/255;
[M1,N1]=size(im_orig);
im=[im_orig(:,1) im_orig im_orig(:,N1)];
im=[im(1,:);im;im(M1,:)];
M2=1000;
N2=1000;
im_new=zeros(M2,N2);

%% Interpolation
for m=1:M2
    for n=1:N2
        mm=fix(m/M2*M1)+1;
        nn=fix(n/N2*N1)+1;
        alpha=m/M2*M1-mm+1;
        beta=n/N2*N1-nn+1;
        im_new(m,n)=(1-alpha)*(1-beta)*im(mm,nn)+alpha*(1-beta)*im(mm+1,nn)...
                +(1-alpha)*beta*im(mm,nn+1)+alpha*beta*im(mm+1,nn+1);
    end
end

%% Result
figure
imagesc(im_orig);
colormap(gray(256));
axis image
figure
imagesc(im_new);
colormap(gray(256));
axis image
