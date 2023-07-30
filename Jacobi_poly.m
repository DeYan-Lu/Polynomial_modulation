function hk = Jacobi_poly(n,a,b)
if n==0 
    syms a b
    hk = 1;
elseif n==1
    syms a b
    hk = [a/2+b/2+1 a/2-b/2];
else    
    hkm2 = zeros(1,n+1); % 當n=0的H(n-2) ; 填滿有n+1個係數
    hkm2(n+1) = 1; %降冪排列
    hkm1 = zeros(1,n+1); %當n=1的H(n-1) ; 填滿有n+1個係數
    hkm1(n) = a/2+b/2+1;  %降冪排列
    hkm1(n+1) = a/2-b/2;
    for k=2:n
        
        hk = zeros(1,n+1);
        for e=n-k+1:2:n  %每隔兩個次方就計算
            hk(e) = ((2*k+a+b-1)*(2*k+a+b)/(2*k*(k+a+b)))*hkm1(e+1)...
                +(((2*k+a+b-1)*(a^2-b^2))/(2*k*(k+a+b)*(2*k+a+b-2)))*hkm1(e)...
                -(((k+a-1)*(k+b-1)*(2*k+a+b))/(k*(k+a+b)*(2*k+a+b-2)))*hkm2(e); %遞迴關係式
        end
        %最後一項為常數項,而H(n-1)有乘以x,所以沒有常數項
       
        hk(n+1) = +(((2*k+a+b-1)*(a^2-b^2))/(2*k*(k+a+b)*(2*k+a+b-2)))*hkm1(n+1)...
                 -(((k+a-1)*(k+b-1)*(2*k+a+b))/(k*(k+a+b)*(2*k+a+b-2)))*hkm2(n+1);
        
        if k<n   %遞迴還沒完時,
            hkm2 = hkm1;
            hkm1 = hk;
        end
        
    end
    
end

