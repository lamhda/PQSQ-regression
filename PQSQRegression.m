function [b, intercept] = PQSQRegression(X, Y, varargin)
%PQSQRegression calculates PQSQ approximation for lasso regression.
%Syntax
%   [b, intercept] = PQSQRegression(X, Y)
%   [b, intercept] = PQSQRegression(X, Y, Name, Value)
%
%Inputs
%   X is numeric matrix with n rows and p columns. Each row represents one
%       observation, and each column represents one predictor (variable). 
%   Y is numeric vector of length n, where n is the number of rows of X.
%       Y(i) is the response to row i of X.
%   Name, Value is one or more pairs of name and value. There are several
%       possible names with corresponding values:
%       'Intervals', intervals serves to specify user defined intervals.
%           intervals is row vector The first element must be zero. By
%           default is created by usage of 'number_of_intervals' and
%           'intshrinkage'. Maximal value M is maximum of absolute value of
%           coefficients for OLS without any penalties multiplied by
%           'intshrinkage'. All other boreders are calcualted as r(i) =
%           M*i^2/p^2, where p is 'number_of_intervals'. 
%       'Number_of_intervals' specifies the number of intervals to
%           automatic interval calculation. Default value is 5. 
%       'intshrinkage', delta serves to specify delta which is coefficient
%           for intervals shrinkage (see argument delta in
%           defineIntervals). Default value is 1 (no shrinkage).
%       'potential' is majorant function for PQSQ. By default it is L1.
%       'Weights' is vector of observation weights.  Must be a vector of
%           non-negative values, of the same length as columns of X.  At
%           least two values must be positive. (default (1/N)*ones(N,1)). 
%
%Return values:
%   b is vector of the fitted coefficients for model. b have dimension Px1,
%       where P = size(X,2) is the number of predictors
%   intercept is intercept of model.

    %Sanity-check of X and Y
    %X must be real valued matrix without Infs and NaNs
    if ~isreal(X) || ~all(isfinite(X(:))) || isscalar(X) || length(size(X))~=2
        error(['Incorrect value for argument "X".'...
            ' It must be real valued matrix without Infs and NaNs']);
    end
    
    %Define number of predictors
    n = size(X, 1);
    
    %Y must be real valued vector without Infs and NaNs and with number of
    %elements equals to n
    if ~isreal(Y) || ~all(isfinite(Y)) || numel(Y)~=n
        error(['Incorrect value for argument "X". It must be real valued',...
            'vector without Infs and NaNs and with number of elements',...
            'equals to number of rows in X']);
    end
    %Transform to column vector.
    Y=Y(:);
        
    %Get optional parameters
    intervals = [];
    numOfInt = 5;
    delta = 1;
    func = @L1;
    weights = [];
    
    %Search Name-Value pairs
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'intervals')
            intervals = varargin{i+1};
        elseif strcmpi(varargin{i},'number_of_intervals')
            numOfInt = varargin{i+1};
        elseif strcmpi(varargin{i},'intshrinkage')
            delta = varargin{i+1};
        elseif strcmpi(varargin{i},'potential')
            func = varargin{i+1};
        elseif strcmpi(varargin{i},'Weights')
            weights = varargin{i+1};
        end
    end
    
    %Sanity-check of parameters and redefine all necessary values

    %weights
    if isempty(weights)
        weights = ones(n,1);
    else
        %Weights must be a vector of nonnegative finite reals with at least two
        %values greater than zero and with number of elements equal to number
        %of rows in X. 
        if ~isreal(weights) || ~isfinite(weights) || sum(weights<0)>0 ||...
                sum(weights>0)<2 || numel(weights)~=n
            error(['Incorrect value for argument "Weights". It must be ',...
                'a vector of nonnegative finite reals with at least two',...
                'values greater than zero and with number of elements equal',...
                ' to number of rows in X.']);
        end
        weights = weights(:);
    end
    %Normalise
    weights = weights / sum(weights);
    
    %Add column vector of ones to the end of X
    X = [X, ones(n, 1)];
    
    
    %Func must be function handler
    if ~isa(func,'function_handle')
        error(['Incorrect value in "potential" argument.'...
            ' It must be function handler']);
    end
    
    %Solve OLS to obtain information for deviations
    A = X;
    A = bsxfun(@times, A, weights);
    b = (A' * X) \ (A' * Y);

    %Calculate deviations for OLS
    d = abs(Y-X*b);
    
    if isempty(intervals)
        %Function has to create intervals by automatic way
        %numOfInt must be positive integer scalar
        if ~isreal(numOfInt) || ~isfinite(numOfInt) || numOfInt < 1
            error(['Incorrect value of "number_of_intervals" argument' ...
                'It must be positive integer scalar']);
        else
            numOfInt = floor(numOfInt);
        end
        %delta has to be positive real scalar
        if ~isreal(delta) || ~isfinite(delta) || delta < 0
            error(['Incorrect value of "intshrinkage" argument' ...
                'It must be positive real scalar']);
        end
            
        pFunc = definePotentialFunction(max(d), numOfInt, func, delta);
    else
        %intervals must contains non nerative values in ascending order.
        %The first value must be zero.
        if intervals(1)~=0 || ~all(isfinite(intervals)) ...
                || any((intervals(2:end)-intervals(1:end-1))<=0)
            error(['Incorrect values in argument intervals: intervals must'...
                ' contains non negative values in ascending order.'...
                ' The first value must be zero.']);
        end
        pFunc.intervals = [intervals, Inf(size(X,2),1)];
        [pFunc.A,pFunc.B] = ...
            computeABcoefficients(intervals, func);
    end
    
    %Main loop
    q = ones(n,1);
    while true
        %Store old values of indeices
        qOld = q;
        %Calculate new deviation absolute values
        d = abs(Y-X*b);
        %Calculate new indices
        q = discretize(d,pFunc.intervals);
        %Check convexity
        if ~any(q-qOld)
            break;
        end
        %Form weights matrix
        d = weights .* (pFunc.A(q))';
        %Calculate new SLAE and solve it
        A = X;
        A = bsxfun(@times, A, d);
        b = (A' * X) \ (A' * Y);
    end
    
    %Form output arguments
    intercept = b(end);
    b = b(1:end-1);
end

function potentialFunction = definePotentialFunction( x,...
    number_of_intervals, potential_function_handle, delta )
%definePotentialFunction defines "uniform in square" intervals for trimming
%threshold x and specified number_of_intervals.
%   x is upper boundary of the interval last but one.
%   number_of_intervals is required number of intervals.
%   potential_function_handle is function handler for coefficients
%       calculation.
%   delta is coefficient of shrinkage which is greater than 0 ang not
%       greater than 1.
%Output argument potentialFunction is structure with three fields:
%   intervals is matrix m-by-number_of_intervals. Each row contains
%       number_of_intervals values of thresholds for intervals and one
%       additional value Inf
%   A and B are the m-by-number_of_intervals matrices with quadratic
%       functions coefficients

    if nargin<4 
        delta = 1;
    end
    
    p=number_of_intervals-1;
    
    %intervals is the product of row and maximal coefficient multiplied by delta:
    intervals = (x * delta) * ((0:p)/p).^2;
    
    potentialFunction.intervals = [intervals, Inf(1)];
    potentialFunction.sqint = potentialFunction.intervals.^2;
    [potentialFunction.A,potentialFunction.B] = ...
        computeABcoefficients(intervals, potential_function_handle);
end

function [A,B] = computeABcoefficients(intervals, potential_function_handle)
%PQSQR_computeABcoefficients calculates the coefficients a and b for
%quadratic fragments of potential function.
%   intervals is the 1-by-K matrix of intervals' boudaries without final
%       infinit boundary.
%   potential_function_handle is a handle of majorant function.

    %Get dimensions of intervals
    p = size(intervals,2);

    %Preallocate memory
    A = zeros(1,p);
    B = zeros(1,p);

    %Calculate value of function all boundaries
    pxk = potential_function_handle(intervals);
    sxk = intervals.^2;

    A(1:p-1) = (pxk(1:p-1)-pxk(2:p))./(sxk(1:p-1)-sxk(2:p));
    B(1:p-1) = (pxk(2:p).*sxk(1:p-1)-pxk(1:p-1).*sxk(2:p))./...
        (sxk(1:p-1)-sxk(2:p));
    B(p) = pxk(p);
end