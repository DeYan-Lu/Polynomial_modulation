function opening=opening(A,B,t)
opening=A;
for i=1:t
    opening=erosion(opening,B);
end
for i=1:t
    opening=dilation(opening,B);
end
end