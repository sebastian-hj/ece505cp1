function [xValues,searchDir,stepLength,errorValues] = ...
    newtonAlgorithm(f,x0,errorValueMax,iterMax)
% Newton
%
%
% J. Sebastián Hurtado J.


% Create output arrays
xValues = zeros([iterMax size(x0,1)]);
searchDir = zeros([iterMax size(x0,1)]);
stepLength = zeros([iterMax 1]);
errorValues = zeros([iterMax 1]);


% Calculate gradient and Hessian matrix
[f,gf,Hf] = gradientHessian(f);

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


% Set initial singularity flag as false
flagSingularity = false;

while errorValue > errorValueMax && iterNumber < iterMax
    
    % Calculate gradient and Hessian using current x
    gfValueOld = gf(xOldCell{:});
    HfValueOld = Hf(xOldCell{:});
    
    % Test for singularity in Hessian matrix
    if( rcond(HfValueOld) < 1e-12 )
        disp('WARNING: Hessian matrix is singular.');
        disp('No optimal value is found.');
        flagSingularity = true;
        break;
    end
    
    % Calculate search direction
    p = HfValueOld\gfValueOld;
    
    % Calculate new function x values
    xNew = xOld - p;
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
    searchDir(iterNumber,:) = gfValueOld';
    stepLength(iterNumber) = 0;
    errorValues(iterNumber) = errorValue;
    
end

% Check for singularity flag
if flagSingularity==true
    
    % Reduce output array size if necessary
    if iterNumber < iterMax
        
        xValues(iterNumber+1:end,:) = '';
        searchDir(iterNumber+1:end,:) = '';
        stepLength(iterNumber+1:end) = '';
        errorValues(iterNumber+1:end) = '';
        
    end
    
else
    
    % Reduce output array size if necessary
    if iterNumber < iterMax
        
        xValues(iterNumber+1:end,:) = '';
        searchDir(iterNumber+1:end,:) = '';
        stepLength(iterNumber+1:end) = '';
        errorValues(iterNumber+1:end) = '';
        
    end
end
    


end