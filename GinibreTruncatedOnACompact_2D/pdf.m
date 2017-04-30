%% Probability density function
% This function provides the probability density function's analytical
% expression needed to perform the rejection sampling and determine the
% points locations.

function pdfiAno = pdf(v, vNorm, i,N, e)

tempsum = 0;
for j = 1:N-i
    tempsum = tempsum + abs(v*conj(e(:,j)))^2;
end

pdfi = 1/i*(vNorm - tempsum);
pdfiAno = matlabFunction(pdfi);
end 