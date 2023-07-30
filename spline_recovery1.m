clear;
%% Image
% filename='C:\Users\DJJ\Documents\Pic\Pepper512c.bmp';
% filename='C:\Users\DJJ\Documents\Pic\Lena512c.bmp';
filename='C:\Users\DJJ\Documents\Pic\Baboon1.bmp';
% filename1='MaxAlpha=2.bmp';
im_orig=double(imread(filename));
im_new=zeros(size(im_orig));
im1=zeros(size(im_orig));
im2=im1;
[M1,N1,ss]=size(im_orig);
im=im_orig; % im=[im_orig im_orig(:,N1)];
nk=2;
ima=zeros(M1,N1+nk*2,ss); ima(:,nk+1:N1+nk,:)=im;
%% Distortion
beta=0.4;
alpha_min=0.2;
m=[1:M1]';
% alpha=alpha_min.^(((M1-m)/(M1-1)).^beta);
% alpha=1+alpha_min-alpha_min.^(((m-1)/(M1-1)).^beta);
alpha=alpha_min+(1-alpha_min)*((m-1)/(M1-1)).^beta;
figure(1)
% subplot(321); plot(alpha)
subplot(321); image(im/256)
nc=round((N1+1)/2);
for m=1:M1
    n21=round((1-nc)*alpha(m));  
    n22=round((N1-nc)*alpha(m));
    n23=[n21:n22];     n24=n23/alpha(m);
    nn=floor(n24);     nn1=nn+nc;
    nx=n24-nn;    
    im_new(m,nc+n23,1)=(1-nx).*ima(m,nn1+nk,1)+nx.*ima(m,nn1+nk+1,1);  
    im_new(m,nc+n23,2)=(1-nx).*ima(m,nn1+nk,2)+nx.*ima(m,nn1+nk+1,2); 
    im_new(m,nc+n23,3)=(1-nx).*ima(m,nn1+nk,3)+nx.*ima(m,nn1+nk+1,3); 
end
subplot(322); image(im_new/255)
%% Recovery
% Linear interpolation
n1=[1-nc:N1-nc];
im_new1=im_new(:,[1:N1,N1],:);
for m=1:M1
    n11=n1*alpha(m);
    n12=floor(n11);     n13=n11-n12; 
    n14=nc+n12;
    im1(m,:,1)=(1-n13).*im_new1(m,n14,1)+n13.*im_new1(m,n14+1,1);
    im1(m,:,2)=(1-n13).*im_new1(m,n14,2)+n13.*im_new1(m,n14+1,2);
    im1(m,:,3)=(1-n13).*im_new1(m,n14,3)+n13.*im_new1(m,n14+1,3);
end
subplot(323); image(im1/255)
Es=sum(sum(sum(im.^2)));
em1=(im1(:,:,1)-im(:,:,1)).^2+(im1(:,:,2)-im(:,:,2)).^2+(im1(:,:,3)-im(:,:,3)).^2;
subplot(324)
image(em1*1); colormap(gray(256))
NRMSE=(sum(sum(em1))/Es)^0.5;
im7=zeros(M1,N1,3);
im_new2=zeros(M1,N1+4,3);
im_new2(:,3:N1+2,:)=im_new;
g1a=1*abs(im_new2(:,3:N1+2,:)-im_new2(:,2:N1+1,:)); %+0*abs(im_new2(:,3:N1+2,:)-im_new2(:,1:N1,:));
g1=g1a(:,:,1)+g1a(:,:,2)+g1a(:,:,3);
g2a=1*abs(im_new2(:,3:N1+2,:)-im_new2(:,4:N1+3,:)); %+0*abs(im_new2(:,3:N1+2,:)-im_new2(:,5:N1+4,:));
g2=g2a(:,:,1)+g2a(:,:,2)+g2a(:,:,3);
gw1=30./(85+1*g1.^0.8);   gw1=gw1(:,[1:N1,N1]);  % gw1=30./(80+g1.^1);
gw2=30./(85+1*g2.^0.8);   gw2=gw2(:,[1:N1,N1]);  % gw2=30./(80+g2.^1);
for m=1:M1
    n11=n1*alpha(m);
    n12=floor(n11);     n13=n11-n12; 
    n14=nc+n12;
    w1=(1-n13).*gw1(m,n14);   w2=n13.*gw2(m,n14+1); 
    w11=w1./(w1+w2);          w21=w2./(w1+w2);    
    im7(m,:,1)=w11.*im_new1(m,n14,1)+w21.*im_new1(m,n14+1,1);
    im7(m,:,2)=w11.*im_new1(m,n14,2)+w21.*im_new1(m,n14+1,2);
    im7(m,:,3)=w11.*im_new1(m,n14,3)+w21.*im_new1(m,n14+1,3);
end
em7=(im7(:,:,1)-im(:,:,1)).^2+(im7(:,:,2)-im(:,:,2)).^2+(im7(:,:,3)-im(:,:,3)).^2;
NRMSE(7)=(sum(sum(em7))/Es)^0.5;
% Spline
im_new2=max(max(im_new(:,:,1),im_new(:,:,2)),im_new(:,:,3));
for m=1:M1
    index1=find(im_new2(m,:),1);
    index2=find(im_new2(m,:),1,'last');
    vn=[index1:index2];
    vn1=(vn-nc)/alpha(m);    
    im2(m,:,1)=spline(vn1,im_new(m,index1:index2,1),n1);
    im2(m,:,2)=spline(vn1,im_new(m,index1:index2,2),n1);
    im2(m,:,3)=spline(vn1,im_new(m,index1:index2,3),n1);
end
im21=im2;
im8=im2*1.1+im1*(-0.1);
im2=min(max(im2,0),255);
im8=min(max(im8,0),255);
em2=(im2(:,:,1)-im(:,:,1)).^2+(im2(:,:,2)-im(:,:,2)).^2+(im2(:,:,3)-im(:,:,3)).^2;
NRMSE(2)=(sum(sum(em2))/Es)^0.5;
em8=(im8(:,:,1)-im(:,:,1)).^2+(im8(:,:,2)-im(:,:,2)).^2+(im8(:,:,3)-im(:,:,3)).^2;
NRMSE(8)=(sum(sum(em8))/Es)^0.5;
% Other Methods
% nearest, pchip, cubic, v5cubic
im3=zeros(M1,N1,3);
im4=zeros(M1,N1,3);
im5=zeros(M1,N1,3);
for m=1:M1
    index1=find(im_new2(m,:),1);
    index2=find(im_new2(m,:),1,'last');
    vn=[index1:index2];
    vn1=(vn-nc)/alpha(m);    
    im3(m,:,1)= interp1(vn1,im_new(m,index1:index2,1),n1,'nearest','extrap');
    im3(m,:,2)= interp1(vn1,im_new(m,index1:index2,2),n1,'nearest','extrap');
    im3(m,:,3)= interp1(vn1,im_new(m,index1:index2,3),n1,'nearest','extrap');
    im4(m,:,1)= interp1(vn1,im_new(m,index1:index2,1),n1,'pchip','extrap');
    im4(m,:,2)= interp1(vn1,im_new(m,index1:index2,2),n1,'pchip','extrap');
    im4(m,:,3)= interp1(vn1,im_new(m,index1:index2,3),n1,'pchip','extrap');
    im5(m,:,1)= interp1(vn1,im_new(m,index1:index2,1),n1,'v5cubic');
    im5(m,:,2)= interp1(vn1,im_new(m,index1:index2,2),n1,'v5cubic');
    im5(m,:,3)= interp1(vn1,im_new(m,index1:index2,3),n1,'v5cubic');
end
im3=min(max(im3,0),255);
im41=im4;
im4=min(max(im4,0),255);
im5=min(max(im5,0),255);
im6=zeros(M1,N1,3);
for m=1:M1
    index1=find(im_new2(m,:),1);
    index2=find(im_new2(m,:),1,'last');
    vn=[index1:index2];
    vn1=(vn-nc)/alpha(m);   
    n11=repmat(n1,[length(vn1),1])-repmat(vn1',[1,N1]);    
    sc=sinc(n11*alpha(m));
    im6(m,:,1)=im_new(m,index1:index2,1)*sc;
    im6(m,:,2)=im_new(m,index1:index2,2)*sc;
    im6(m,:,3)=im_new(m,index1:index2,3)*sc;
end
im6=min(max(im6,0),255);
im9=im21*0.7+im41*(0.3);
im9=min(max(im9,0),255);
gw=0.6+tanh((g1+g2-30)*0.1)*0.2;
im10=zeros(M1,N1,3);
im10(:,:,1)=im21(:,:,1).*gw+im41(:,:,1).*(1-gw);
im10(:,:,2)=im21(:,:,2).*gw+im41(:,:,2).*(1-gw);
im10(:,:,3)=im21(:,:,3).*gw+im41(:,:,3).*(1-gw);
im10=min(max(im10,0),255);
gw1=1.25+tanh((g1+g2-15)*0.2)*0.35;
im11=zeros(M1,N1,3);
im11(:,:,1)=im21(:,:,1).*gw1+im1(:,:,1).*(1-gw1);
im11(:,:,2)=im21(:,:,2).*gw1+im1(:,:,2).*(1-gw1);
im11(:,:,3)=im21(:,:,3).*gw1+im1(:,:,3).*(1-gw1);
im11=min(max(im11,0),255);
subplot(325); image(im2/255)
subplot(326)
image(em2*1); colormap(gray(256))
em3=(im3(:,:,1)-im(:,:,1)).^2+(im3(:,:,2)-im(:,:,2)).^2+(im3(:,:,3)-im(:,:,3)).^2;
NRMSE(3)=(sum(sum(em3))/Es)^0.5;
em4=(im4(:,:,1)-im(:,:,1)).^2+(im4(:,:,2)-im(:,:,2)).^2+(im4(:,:,3)-im(:,:,3)).^2;
NRMSE(4)=(sum(sum(em4))/Es)^0.5;
em5=(im5(:,:,1)-im(:,:,1)).^2+(im5(:,:,2)-im(:,:,2)).^2+(im5(:,:,3)-im(:,:,3)).^2;
NRMSE(5)=(sum(sum(em5))/Es)^0.5;
em6=(im6(:,:,1)-im(:,:,1)).^2+(im6(:,:,2)-im(:,:,2)).^2+(im6(:,:,3)-im(:,:,3)).^2;
NRMSE(6)=(sum(sum(em6))/Es)^0.5;
em9=(im9(:,:,1)-im(:,:,1)).^2+(im9(:,:,2)-im(:,:,2)).^2+(im9(:,:,3)-im(:,:,3)).^2;
NRMSE(9)=(sum(sum(em9))/Es)^0.5;
em10=(im10(:,:,1)-im(:,:,1)).^2+(im10(:,:,2)-im(:,:,2)).^2+(im10(:,:,3)-im(:,:,3)).^2;
NRMSE(10)=(sum(sum(em10))/Es)^0.5;
em11=(im11(:,:,1)-im(:,:,1)).^2+(im11(:,:,2)-im(:,:,2)).^2+(im11(:,:,3)-im(:,:,3)).^2;
NRMSE(11)=(sum(sum(em11))/Es)^0.5;
NRMSE*10
figure(2)
subplot(321); image(im8/255)
subplot(322)
image(em8*1); colormap(gray(256))
subplot(323); image(im9/255)
subplot(324)
image(em9*1); colormap(gray(256))
subplot(325); image(im11/255)
subplot(326)
image(em11*1); colormap(gray(256))

