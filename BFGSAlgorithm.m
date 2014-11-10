function [xValues,searchDir,stepLength,errorValues] = ...
    BFGSAlgorithm(f,x0,errorValueMax,iterMax)
% BFGS
%
%
% J. Sebastián Hurtado J.


% Create output arrays
xValues = zeros([iterMax size(x0,1)]);
searchDir = zeros([iterMax size(x0,1)]);
stepLength = zeros([iterMax 1]);
errorValues = zeros([iterMax 1]);


% Calculate gradient
[f,gf] = gradientHessian(f);

% Set inital value as old value
xOld = x0;
% Convert to cell
xOldCell = num2cell(xOld);

% Initial error and number of iteration
errorValue = 1;
iterNumber = 1;


% Add inital values to output array
xValues(iterNumber,:) = xOld';
errorValues(iterNumber) = norm(gf(xOldCell{:}))/(1 + abs(f(xOldCell{:})));


% Calculate inital B matrix
BOld = eye(size(xOld,1));


while errorValue > errorValueMax && iterNumber < iterMax
    
    % Calculate gradient
    gfValueOld = gf(xOldCell{:});
    
    % Find p vector
    p = BOld\(-gfValueOld);
    
    % Find optimum alpha
    %alphaValue = simpleLineSearch(f,p,xOld);
    alphaValue = simpleLineSearch(f,p,gfValueOld,xOld,0.3);
    
    % Calculate new function x values
    xNew = xOld + alphaValue*p;
    % Convert to cell
    xNewCell = num2cell(xNew);
    
    % Calculate 's' and 'y' vectors
    s = xNew - xOld;
    y = gf(xNewCell{:}) - gf(xOldCell{:});
    
    % Hessian approximation
    BNew = BOld - ( (BOld*s)*(BOld*s)' )/( s'*BOld*s ) + (y*y')/(y'*s);
    
    % Calculate error
    fValueNew = f(xNewCell{:});
    gfValueNew = gf(xNewCell{:});
    errorValue = norm(gfValueNew)/(1 + abs(fValueNew));
    
    % Update value
    xOld = xNew;
    % Convert to cell
    xOldCell = num2cell(xOld);
    
    % Update Hessian approximation
    BOld = BNew;
    
    % Update iteration
    iterNumber = iterNumber + 1;
    
    % Add to output array
    xValues(iterNumber,:) = xNew';
    searchDir(iterNumber,:) = p';
    stepLength(iterNumber) = alphaValue;
    errorValues(iterNumber) = errorValue;
    
end

% Reduce output array size if necessary
if iterNumber < iterMax
   
    xValues(iterNumber+1:end,:) = '';
    searchDir(iterNumber+1:end,:) = '';
    stepLength(iterNumber+1:end) = '';
    errorValues(iterNumber+1:end) = '';
    
end

end
