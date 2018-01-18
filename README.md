# PQSQ-regression
The method of PQSQ regression is presented. This method allow us quickly and with high accuracy fit regresion models other than L2 regression.

Theory of PQSQ approximation can be found in [wiki](https://github.com/Mirkes/PQSQ-regression/wiki)

Current MatLab implementation includes several files:

L1 is majorant function for L_1 norm.

L1_5 is majorant function for L_1 + L_2 norm.

L2 is majorant function for L_2 norm.

LLog is natural logarithm majorant function.

LSqrt is square root or L_0.5 quasi norm majorant function.

PQSQRegression fits regression model with PQSQ penalty of absolute deviations.

<pre>
Syntax
   [b, intercept] = PQSQRegression(X, Y)
   [b, intercept] = PQSQRegression(X, Y, Name, Value)


Example

 %Decomment nex line for reproducibility
 %rand('twister',2323); randn('state',1865);
 
 %Create test data
 %True values of coefficients
 BTrue = [1, -2];
 nDots = 100;
 X = rand(createNew, 1);
 Y = BTrue(1) + X * BTrue(2) + 0.05 * randn(nDots, 1);
 %Add some outliers by transfer nOutlier points by vector [-1, -1]
 nOutlier = floor(nDots / 20);
 X(1:nOutlier) = X(1:nOutlier) - 1;
 Y(1:nOutlier) = Y(1:nOutlier) - 1;
 %Array of axis limits
 lims = [-1, 1, -2, 1.1];
 %Prepare plot
 figure;
 %Draw data points
 plot(X, Y, 'bs');
 hold on;
 %Draw true regression line
 drawPlot(BTrue, lims, 'r-');
 
 %Calculate and draw L2 regression
 BS = [ones(length(X), 1), X] \ Y;
 drawPlot(BS, lims, 'k-');
 
 %PQSQ L2
 BP2 = BS;
 [BP2(2), BP2(1)] = PQSQRegression(X, Y, 'potential', @L2, 'intshrinkage',... 
         0.5, 'Number_of_intervals', 2);
 drawPlot(BP2, lims, 'b-');
 
 %PQSQ L1
 BP1 = BS;
 [BP1(2), BP1(1)] = PQSQRegression(X, Y);
 drawPlot(BP1, lims, 'c-');
 
 %PQSQ L0.5
 BP05 = BS;
 [BP05(2), BP05(1)] = PQSQRegression(X, Y, 'potential', @LSqrt);
 drawPlot(BP05, lims, 'm-');
 
 axis(lims);
 grid on 
 legend('Data', 'Pure', 'L-2', 'PQSQ L-2', 'PQSQ L-1', 'PQSQ L-0.5');
 
 
Inputs
   X is numeric matrix with n rows and p columns. Each row represents one
       observation, and each column represents one predictor (variable). 
   Y is numeric vector of length n, where n is the number of rows of X.
       Y(i) is the response to row i of X.
   Name, Value is one or more pairs of name and value. There are several
       possible names with corresponding values:
       'Intervals', intervals serves to specify user defined intervals.
           intervals is row vector The first element must be zero. By
           default is created by usage of 'number_of_intervals' and
           'intshrinkage'. Maximal value M is maximum of absolute value of
           coefficients for OLS without any penalties multiplied by
           'intshrinkage'. All other boreders are calcualted as r(i) =
           M*i^2/p^2, where p is 'number_of_intervals'. 
       'Number_of_intervals' specifies the number of intervals to
           automatic interval calculation. Default value is 5. 
       'intshrinkage', delta serves to specify delta which is coefficient
           for intervals shrinkage (see argument delta in
           defineIntervals). Default value is 1 (no shrinkage).
       'potential' is majorant function for PQSQ. By default it is L1.
       'Weights' is vector of observation weights.  Must be a vector of
           non-negative values, of the same length as columns of X.  At
           least two values must be positive. (default (1/N)*ones(N,1)). 

Return values:
   b is vector of the fitted coefficients for model. b have dimension Px1,
       where P = size(X,2) is the number of predictors
   intercept is intercept of model.
</pre>