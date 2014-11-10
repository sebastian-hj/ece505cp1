function [fFun,gfFun,HfFun] = gradientHessian(f)
% Calculate gradient and Hessian
%
% Input:
%   - f: function in char class, e.g., 'x1 + x2^2'
%
% J. Sebastián Hurtado J.

% Convert string to symbolic
fSym = sym(f);
% Find symbolic variables in expression
fVars = symvar(fSym);

% Create array for gradient
gf = sym(zeros(size(fVars,2),1));
% Calculate gradient
for k = 1:size(fVars,2)
    gf(k,1) = diff(fSym,fVars(k));
end

% Calculate Hessian
Hf = jacobian(gf,fVars);

% Convert to anonymous function using variables from function
fFun = matlabFunction(fSym,'vars',fVars);
gfFun = matlabFunction(gf,'vars',fVars);
HfFun = matlabFunction(Hf,'vars',fVars);