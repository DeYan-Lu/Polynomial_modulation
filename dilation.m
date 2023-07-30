function dil=dilation(A,B)
[M,N]=size(A);
dil=zeros(size(A));
for i=1:M
    for j=1:N
        if A(i,j)==1
            for k=1:size(B,1)
                if (i+B(k,1))>0 && (i+B(k,1))<=M && (j+B(k,2))>0 && (j+B(k,2))<=N
                    dil(i+B(k,1),j+B(k,2))=1;
                end
            end
        end
    end
end
end