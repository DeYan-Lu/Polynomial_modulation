function P = binomial(n,a,b)
%------------------------------------------------------
%
% function P = binomial(n,a,b)
%
% Returns the polynomial coefficients p(i) for
%
% P(x) = (ax + b)^n
%      = p(1)x^n + p(2)x^(n-1) + ... + p(n)x + p(n+1)
%
%------------------------------------------------------
% First binomial coefficients for (x+1)^n
if n==0, P=1; return, end
p = zeros(n,n+1);
p(1,1:2) = [1 1];
for j = 2:n
    p(j,1) = 1;
    for k = 2:n
        p(j,k) = p(j-1,k-1) + p(j-1,k);
    end
    p(j,n+1) = 1;
end
% Adjust for a and b
for k = 1:n+1;
    P(k) = p(n,k)*a^(n-k+1)*b^(k-1);
end
end


