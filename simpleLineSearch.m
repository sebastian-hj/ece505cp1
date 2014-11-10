function alphaValue = simpleLineSearch(f,p,gf,x,u)
% Simple line search method
%
%
% J. Sebastián Hurtado J.

% Start with alphaValue = 1;
alphaValue = 1;

% Convert to cell
xCell = num2cell(x);

% Alpha cell
xAlphaCell = num2cell(x + alphaValue*p);

% Calculate initial value for function
fAlpha = f(xAlphaCell{:});

% Computa Armijo condition
fArmijo = f(xCell{:}) + u*alphaValue*p'*gf;

while fAlpha > fArmijo
    
    % Reduce alphaValue by half
    alphaValue = alphaValue/2;
    
    % Compute alpha
    xAlphaCell = num2cell(x + alphaValue*p);
    fAlpha = f(xAlphaCell{:});
    
    % Compute Armijo
    fArmijo = f(xCell{:}) + u*alphaValue*p'*gf;
    
end


end