%% This function allows to compute versor e, needed in the calculation of the pdf.

function ecol = computeE(v, N, i, X, e)

vAno = matlabFunction(v);
tempsum = 0;

for j = 1:N-i
    tempsum = tempsum + (conj(e(:, j))'.*feval(vAno, X(i)))*e(:, j);
end
wi = feval(vAno, X(i)) - tempsum;

ecol = (wi/norm(wi))';

end