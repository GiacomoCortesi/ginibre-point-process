
function pdfiAno = pdf(v, vNormSquare, i, N, e)

tempsum = 0;
for j = 1:N-i
    tempsum = tempsum + abs(v*conj(e(:,j)))^2;
end

pdfi = 1/i*(vNormSquare - tempsum);
pdfiAno = matlabFunction(pdfi);


%% Implementation notes for Gaussian distribution (not needed)
% out = N+a+x;
% Now u have to generate random numbers given a specific probability density function
% The case XN should reduce to a gaussian distribution (to be proved):
% X = 1:N;
% MU = 2;
% SIGMA = 5;
% X(N) = normrnd(MU, SIGMA);
% x = -20:.1:20;
% normal = normpdf(x, MU, SIGMA);
% plot(x, normal);
% pdfN = norm(v)
% f = @(q) 1/20;
% X = sampleDist(@(x) 1/sqrt(2*pi) *exp(-x.^2/2),1/sqrt(2*pi),1e6,[-5,5],true);
end 