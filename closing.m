function closing=closing(A,B,t)
closing=A;
for i=1:t
    closing=dilation(closing,B);
end
for i=1:t
    closing=erosion(closing,B);
end
end