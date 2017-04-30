%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                                      %%%%%%%
%%%%%%%        MONO-DIMENSIONAL GINIBRE POINT PROCESS        %%%%%%%
%%%%%%%                                                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation of a Ginibre Point Process on the mono-dimensional plane.
% The user sets:
% - N: number of nodes
% - a: radius of the ball


%% Start fresh with the simulation %%
% Clear variables and functions from memory.
clear; 
% Close all the open figure windows
close all; 
% Clear command window
clc;

%% Let the user insert the value of a and N from the command window %%
    a = input('Enter the value of a: ');
    N = input('Enter the value of N: ');

%% Algorithm Implementation %%

% Note1: When a variable contains the word 'Ano' it means that it's an
% anonimouos function. In fact Matlab has a bit tricky way of dealing with
% functions, if more information are needed refer to: Symbolic and
% anonimouos functions in the Matlab documentation.
%
% Note2: The Rejection Sampling algorithm has been taken from the web. It
% has been tried and it worked fine for the Gaussian and Uniform pdf cases.
%
% Note3: The second input parameter of the acceptance-rejection algorithm,
% that is M, should be the maximum value that can be assumed by the pdf. In
% the code is set to one but it has been evaluated intuitevely with no
% analytical proof. A possible improvement could be verifying the maximum
% numerically through Matlab.
% 

% fi function defined using Symbolic variables.
syms t;
syms z;
syms k;
fi = symfun((N/(pi*a^2*int(exp(-t)*t^(k+1-1), 0, N)))*(exp(-(N*(abs(z)^2))/(2*a^2)))*(((N*z)/a^2)^k),[z,k]);

% This compute the vector v of size N, each element is a fi function.
% In order to reduce the computational time, v should be
% initialized before the for cycle. However, I couldn't figure out how to do
% it since it will be filled with Symbolic functions.
for k = 0:N-1
    temp = fi(z, k);
    v(k+1) = temp;
end

% This is in order to compute the square value of the norm of vector v.
vNormSquare = 0;
for i = 1:N
     vNormSquare = vNormSquare + v(i)^2;
end


pdfN = vNormSquare/N;
pdfNAno = matlabFunction(vNormSquare);

%% pdfN normalisation (so that the integral gives one)
pdfNNorm = pdfN/int(pdfN,-a,a);
pdfNNormAno = matlabFunction(pdfNNorm);

%% Assignment of the first sample generated
Xvec = sampleDist(pdfNAno,10,100000,[-a,a]);
XN = Xvec(1);
vXn = 1:N;

for i = 1:N
    vAno = matlabFunction(v(i));
    vXn(i) = vAno(XN);
end

vNormAno = matlabFunction(vNormSquare);
vNormXn = vNormAno(XN);

e1 = vXn / sqrt(vNormXn);
% Check: e1 is a versor (vector with norm 1)
norm(e1);

e = zeros(N, N);
X(N) =  XN;
e(:,1) = e1;

for i=N-1:-1:1
    pdfiAno = pdf(v, vNormSquare, i, N, e);
    pdfi = sym(pdfiAno);
    
    %% pdfi normalisation (so that the integral gives one)
    pdfiNorm = pdfi/int(pdfi,-a,a);
    pdfiNormAno = matlabFunction(pdfiNorm);
    
    Xvec = sampleDist(pdfiAno,10,100000,[-a,a]);
    X(i) = Xvec(1);
    ecol = computeE(v, N, i, X, e);
    e(:,N-i+1) = ecol;
end

% Vector X contains now the samples representing the point locations along
% a line

plot(X, zeros(length(X)), 'ok');
title('Mono-dimensional Ginibre Truncated on a Compact')

