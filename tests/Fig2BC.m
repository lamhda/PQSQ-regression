% Figure 2B
% Data
points = [100, 1000, 10000, 100000];
timeL1 = [0.7152, 0.5709, 22.9928, 475.8805];
timePL2 = [0.1931, 0.0425, 1.0700, 8.8527];
timePL1 = [0.2778, 0.0524, 1.9054, 16.9613];
timePL05 = [0.3059, 0.0598, 1.6563, 16.7613];
angleL1 = [0.1896, 0.5761, 0.4811, 0.4606];
anglePL2 = [0.3320, 0.0839, 0.0120, 0.0074];
anglePL1 = [0.3320, 0.0906, 0.0165, 0.0080];
anglePL05 = [0.3320, 0.0917, 0.0172, 0.0081];

width = 2;
%Draw figure 2B
figure;
loglog(points, timeL1, 'k-', 'LineWidth', width);
hold on;
loglog(points, timePL2, 'r--', 'LineWidth', width);
loglog(points, timePL1, 'b:', 'LineWidth', width);
loglog(points, timePL05, 'g-.', 'LineWidth', width);
legend('L-1', 'PQSQ L-2', 'PQSQ L-1', 'PQSQ L-0.5','Location','northwest');
set(gca,'ytick',[0.1, 10, 1000]);
xlabel('Number of points');
set(gca,'fontsize', 14);
title('Time spent, in seconds','fontsize', 16);
set(gcf,'pos',[500,10,560, 280]);

%Draw figure 2C
figure;
loglog(points, angleL1, 'k-', 'LineWidth', width);
hold on;
loglog(points, anglePL2, 'r--', 'LineWidth', width);
loglog(points, anglePL1, 'b:', 'LineWidth', width);
loglog(points, anglePL05, 'g-.', 'LineWidth', width);
legend('L-1', 'PQSQ L-2', 'PQSQ L-1', 'PQSQ L-0.5','Location','southwest');
xlabel('Number of points');
set(gca,'fontsize', 14);
title('Angular error, in degrees','fontsize', 16);
set(gcf,'pos',[500,10,560, 280]);