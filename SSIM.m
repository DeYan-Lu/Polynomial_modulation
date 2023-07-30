function S = SSIM(X,Y,c1,c2)
% X: target, Y: input for comparison
[M,N]=size(X);
if M~=1
    L=255;
    ax=mean(mean(X));
    ay=mean(mean(Y));
    vx2=sum(sum(abs(X-ax).^2))/M/N;
    vy2=sum(sum(abs(Y-ay).^2))/M/N;
    vxy=sum(sum((X-ax).*(Y-ay)))/M/N;
    S=(2*ax*ay+c1*L)/(ax^2+ay^2+c1*L)*(2*vxy+c2*L)/(vx2+vy2+c2*L);
else
    L=max(X);
    ax=mean(X);
    ay=mean(Y);
    vx2=sum(abs(X-ax).^2)/N;
    vy2=sum(abs(Y-ay).^2)/N;
    vxy=sum((X-ax).*(Y-ay))/N;
    S=(2*ax*ay+c1*L)/(ax^2+ay^2+c1*L)*(2*vxy+c2*L)/(vx2+vy2+c2*L);
end
