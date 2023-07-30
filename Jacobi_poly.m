function hk = Jacobi_poly(n,a,b)
if n==0 
    syms a b
    hk = 1;
elseif n==1
    syms a b
    hk = [a/2+b/2+1 a/2-b/2];
else    
    hkm2 = zeros(1,n+1); % ��n=0��H(n-2) ; �񺡦�n+1�ӫY��
    hkm2(n+1) = 1; %�����ƦC
    hkm1 = zeros(1,n+1); %��n=1��H(n-1) ; �񺡦�n+1�ӫY��
    hkm1(n) = a/2+b/2+1;  %�����ƦC
    hkm1(n+1) = a/2-b/2;
    for k=2:n
        
        hk = zeros(1,n+1);
        for e=n-k+1:2:n  %�C�j��Ӧ���N�p��
            hk(e) = ((2*k+a+b-1)*(2*k+a+b)/(2*k*(k+a+b)))*hkm1(e+1)...
                +(((2*k+a+b-1)*(a^2-b^2))/(2*k*(k+a+b)*(2*k+a+b-2)))*hkm1(e)...
                -(((k+a-1)*(k+b-1)*(2*k+a+b))/(k*(k+a+b)*(2*k+a+b-2)))*hkm2(e); %���j���Y��
        end
        %�̫�@�����`�ƶ�,��H(n-1)�����Hx,�ҥH�S���`�ƶ�
       
        hk(n+1) = +(((2*k+a+b-1)*(a^2-b^2))/(2*k*(k+a+b)*(2*k+a+b-2)))*hkm1(n+1)...
                 -(((k+a-1)*(k+b-1)*(2*k+a+b))/(k*(k+a+b)*(2*k+a+b-2)))*hkm2(n+1);
        
        if k<n   %���j�٨S����,
            hkm2 = hkm1;
            hkm1 = hk;
        end
        
    end
    
end

