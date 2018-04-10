%Decomment nex line for reproducibility
%rand('twister',2323); randn('state',1865);
%Create test data or load previously created
if createNew > 0
    %True values of coefficients
    BTrue = [1, -2]; 
    X = rand(createNew, 1);
    Y = BTrue(1) + X * BTrue(2) + 0.05 * randn(createNew, 1);
    %Add some outliers by transfer nOutlier points by vector [-1, -1]
    X(1:nOutlier) = X(1:nOutlier) - 1;
    Y(1:nOutlier) = Y(1:nOutlier) - 1;
else
    load('AtrificialTestSet100.mat');  %    100 points     5 outliers
    %load('AtrificialTestSet1k.mat');   %  1,000 points    50 outliers
    %load('AtrificialTestSet10k.mat');  % 10,000 points   500 outliers
    %load('AtrificialTestSet100k.mat');  %100,000 points 5,000 outliers
end
%Array of axis limits
lims = [-1, 1, -2, 1.1];
% Number of repetitions for time estimations
%It is good idea to use at least 
%10,000 for     100 points. 
%10,000 for   1,000 points
% 1,000 for  10,000 points
%   100 for 100,000 points
nReps = 1000; 
times = zeros(5,1);

%Prepare plot
figure;
%Draw data points
plot(X, Y, 'bs');
hold on;
%Draw true regression line
drawPlot(BTrue, lims, 'r-');

%Calculate and draw L2 regression
tic;
for k = 1:nReps
    BS = [ones(length(X), 1), X] \ Y;
end
timeTest(1) = toc;
drawPlot(BS, lims, 'k--');

%Calculate L1 regression and draw it
tic;
for k = 1:nReps
    BL1 = L1LinearRegression(X,Y); 
end
timeTest(2) = toc;
drawPlot(BL1, lims, 'k:');

%PQSQ L2
BP2 = BS;
tic;
for k = 1:nReps
    [BP2(2), BP2(1)] = PQSQRegression(X, Y, 'potential', @L2, 'intshrinkage',... 
        0.5, 'Number_of_intervals', 2);
end
timeTest(3) = toc;
drawPlot(BP2, lims, 'b-.');

%PQSQ L1
BP1 = BS;
tic;
for k = 1:nReps
    [BP1(2), BP1(1)] = PQSQRegression(X, Y);
end
timeTest(4) = toc;
drawPlot(BP1, lims, 'c-', 2);

%PQSQ L0.5
BP05 = BS;
tic;
for k = 1:nReps
    [BP05(2), BP05(1)] = PQSQRegression(X, Y, 'potential', @LSqrt);
end
timeTest(5) = toc;
drawPlot(BP05, lims, 'm--', 2);

axis(lims);
grid on 
legend('Data', 'Pure', 'L-2', 'L-1', 'PQSQ L-2', 'PQSQ L-1', 'PQSQ L-0.5');

BB = [BS'; BL1'; BP2'; BP1'; BP05'];

set(gca,'fontsize', 14);
set(gcf,'pos',[10,10,560, 560]);