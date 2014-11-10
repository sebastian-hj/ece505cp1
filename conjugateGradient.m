function [xValues,searchDir,stepLength,errorValues] = ...
    conjugateGradient(f,x0,errorValueMax,iterMax)
% Conjugate Gradient
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

% Initial d
dVector = -gf(xOldCell{:});

% Add inital values to output array
xValues(iterNumber,:) = xOld';
errorValues(iterNumber) = norm(gf(xOldCell{:}))/(1 + abs(f(xOldCell{:})));
%searchDir(iterNumber,:) = dVector';

while errorValue > errorValueMax && iterNumber < iterMax
    
    % Calculate gradient and Hessian using current x
    gfValueOld = gf(xOldCell{:});
    HfValueOld = Hf(xOldCell{:});
    
    % Calculate alpha value
    alphaValue = - ( dVector'*gfValueOld )/( dVector'*HfValueOld*dVector );
    
    % Calculate new x values
    xNew = xOld + alphaValue*dVector;
    % Convert to cell
    xNewCell = num2cell(xNew);
    
    % Calculate error
    fValueNew = f(xNewCell{:});
    gfValueNew = gf(xNewCell{:});
    errorValue = norm(gfValueNew)/(1 + abs(fValueNew));
    
    % Calculate beta value
    beta = ( dVector'*HfValueOld*gfValueNew )/( dVector'*HfValueOld*dVector );
    
%     % Calculate new d value
%     dVector = -gfValueNew + beta*dVector;
    
    % Update value
    xOld = xNew;
    % Convert to cell
    xOldCell = num2cell(xOld);
    
    % Update iteration
    iterNumber = iterNumber + 1;
    
    % Add to output array
    xValues(iterNumber,:) = xNew';
    searchDir(iterNumber,:) = dVector';
    stepLength(iterNumber) = alphaValue;
    errorValues(iterNumber) = errorValue;
    
	% Calculate new d value
    dVector = -gfValueNew + beta*dVector;
    
end

% Reduce output array size if necessary
if iterNumber < iterMax
    
    xValues(iterNumber+1:end,:) = '';
    searchDir(iterNumber+1:end,:) = '';
    stepLength(iterNumber+1:end) = '';
    errorValues(iterNumber+1:end) = '';
    
end


end