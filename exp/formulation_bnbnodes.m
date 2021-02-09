clc; clear;

data = dlmread('formulation_bnbnodes.txt');

figfile_all = 'formulation_bnbnode_all.eps';
figfile_best = 'formulation_bnbnode_best.eps';

[out, idx] = sort(data, 2);

N = length(data(:,1));

% ... Y-axis (one step per problem)
Y = [1/N:1/N:1];

% ... sort data for different methods
ours = sort ( min(data(:,1:4), [], 2) );
cozads = sort ( min(data(:,5:8), [], 2) ) ;

% ... plot the stairs (change axis below to zoom in/out)
fig=figure(1); clf;
set(gca,'FontSize',15);
axis([0 150000 0.5 1 ]); hold on;
stairs(ours,Y,'k-','LineWidth',2);  %set(gca, 'YScale', 'log');
stairs(cozads,Y,'k--','LineWidth',2); 
set(gca,'XTick',0:50000:150000);
xlabel('BnBnodes');
ylabel('Fraction of Problems');
legend('Best of Imp','Best of Coz','Location','SouthEast');
saveas(fig, 'formulation_bnbnode_best.eps');
saveas(fig, 'formulation_bnbnode_best.png');

fig=figure(2); clf;
set(gca,'FontSize',15);
axis([0 150000 0.5 1 ]); hold on;
lindwidth = 1.5;
markersize = 5;
stairs(sort(data(:,1)),Y,'k-o','LineWidth', lindwidth, 'MarkerSize', 4);
stairs(sort(data(:,2)),Y,'k-x','LineWidth', lindwidth, 'MarkerSize', markersize);
stairs(sort(data(:,3)),Y,'k-d','LineWidth', lindwidth, 'MarkerSize', markersize);
stairs(sort(data(:,4)),Y,'k-h','LineWidth', lindwidth, 'MarkerSize', markersize);
stairs(sort(data(:,5)),Y,'k:o','LineWidth', lindwidth, 'MarkerSize', 4);
stairs(sort(data(:,6)),Y,'k:x','LineWidth', lindwidth, 'MarkerSize', markersize);
stairs(sort(data(:,7)),Y,'k:d','LineWidth', lindwidth, 'MarkerSize', markersize);
stairs(sort(data(:,8)),Y,'k:h','LineWidth', lindwidth, 'MarkerSize', markersize);
set(gca,'XTick',0:50000:150000);
xlabel('BnBnodes');
ylabel('Fraction of Problems');
legend('Imp-F','Imp-R','Imp-S','Imp-N','Coz-F','Coz-R','Coz-S','Coz-N','Location','NorthEast');
saveas(fig, 'formulation_bnbnode_all.eps')
saveas(fig, 'formulation_bnbnode_all.png')
