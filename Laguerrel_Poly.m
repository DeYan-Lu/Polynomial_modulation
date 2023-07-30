function hk = Laguerrel_Poly(n,a)
if n==0 
    hk = 1;
elseif n==1 
    hk = [-1 a+1];
else    
    hkm2 = zeros(1,n+1); % ��n=0��H(n-2) ; �񺡦�n+1�ӫY��
    hkm2(n+1) = 1; %�����ƦC
    hkm1 = zeros(1,n+1); %��n=1��H(n-1) ; �񺡦�n+1�ӫY��
    hkm1(n) = -1;  %�����ƦC
    hkm1(n+1) = a+1;
    for k=2:n      
        hk = zeros(1,n+1);
        for e=n-k+1:1:n  %�C�j��Ӧ���N�p��
            hk(e) = (-1/k)*hkm1(e+1)...
                +((2*k-1+a)/k)*hkm1(e)...
                -((k-1+a)/k)*hkm2(e); %���j���Y��
        end
        %�̫�@�����`�ƶ�,��H(n-1)�����Hx,�ҥH�S���`�ƶ�
       
        hk(n+1) = +((2*k-1+a)/k)*hkm1(n+1)...
                  -((k-1+a)/k)*hkm2(n+1);       
        if k<n   %���j�٨S����,
            hkm2 = hkm1;
            hkm1 = hk;
        end
        
    end
end
end

