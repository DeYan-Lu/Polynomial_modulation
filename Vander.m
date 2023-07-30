 function a = Vander(x,y,N)
 %
 % x  column vector of domain data points
 % y  column vector of range data points
 %
 % a  column vector of coefficients for the
 %    interpolating polynomial p
 %
 n = length(x);
 V = ones(n,N+1);
 for j=2:N+1
   %Set up column j
   V(:,j) = x'.*V(:,j-1);
 end
 %Use the built in
 a= mldivide(V,y') ;
 a=fliplr(a');

 end
 