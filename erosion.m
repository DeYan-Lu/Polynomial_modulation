function ero=erosion(A,B)
% Use Erosion Dilation Duality
ero=1-(dilation(1-A,-B));
end