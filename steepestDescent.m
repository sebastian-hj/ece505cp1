function [xValues,searchDir,stepLength,errorValues] = ...
    steepestDescent(f,x0,errorValueMax,iterMax)
% Steepest descent method
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

while errorValue > errorValueMax && iterNumber < iterMax
   
    % Calculate gradient value
    gfValueOld = gf(xOldCell{:});
    
    % Find optimum alpha
    %alphaValue = simpleLineSearch(f,-gfValueOld,xOld);
    alphaValue = simpleLineSearch(f,-gfValueOld,gfValueOld,xOld,0.3);
    
    % Calculate new function x values
    xNew = xOld - alphaValue*gfValueOld;
    % Convert to cell
    xNewCell = num2cell(xNew);
    
    % Calculate error
    fValueNew = f(xNewCell{:});
    gfValueNew = gf(xNewCell{:});
    errorValue = norm(gfValueNew)/(1 + abs(fValueNew));
    
    % Update value
    xOld = xNew;
    % Convert to cell
    xOldCell = num2cell(xOld);
    
    % Update iteration
    iterNumber = iterNumber + 1;
    
    % Add to output array
    xValues(iterNumber,:) = xNew';
    searchDir(iterNumber,:) = -gfValueOld';
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