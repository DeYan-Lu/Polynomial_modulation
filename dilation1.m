function dil=dilation1(A,B,center)
[M,N]=size(A);
dil=zeros(size(A));
[fx,fy]=find(B);
fx=fx-center(1);  fy=fy-center(2);  
L=length(fx);
for c=1:L
    if fx(c)<0
        vx=[ones(1,abs(fx(c))),1:M+fx(c)];
    else
        vx=[1+fx(c):M,M*ones(1,fx(c))];
    end
    if fy(c)<0
        vy=[ones(1,abs(fy(c))),1:N+fy(c)];
    else
        vy=[1+fy(c):N,N*ones(1,fy(c))];
    end
    dil=max(dil,A(vx,vy));
end

end