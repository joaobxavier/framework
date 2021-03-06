% scDrawFigure6

% get the simulation data for sim #3
addpath('../ModelResultBrowser');
basedir = 'E:\jxavier';
resultPath = [basedir '/results'];
dirName = [resultPath '/AnaerobicFeedDO40Precision0_01MaxF0_95_CN50'];

% import the experimental data 
experimental = xlsread([basedir...
    '\granular sludge\with-coathors-comments\after-NH4-16052006\figures\figure6\experimental-data-only.xls']);
% experimental series

timeExp = experimental(:, 1) / 60;
nh4Exp = experimental(:, 2) / 1000;
no2Exp = experimental(:, 3) / 1000;
no3Exp = experimental(:, 4) / 1000;
po4Exp = experimental(:, 5) / 1000;
acExp = experimental(:, 6) / 1000;


% defining the colors
substrateColor = [0 0 1];
oxygenColor = [0 0 0];
ammoniumColor = [0 1 0];
nitriteColor = [1 1 0];
nitrateColor = [0 0 0];
phosphateColor = [1 0 0];

activePAOColor = phosphateColor;
phbPAOColor = [0 0 1];
polypPAOColor = [1 0 0];
glycogenPAOColor = [0 1 1];
totalInertColor = [0.5 0.5 0.5];
activeNHColor = ammoniumColor;
activeNOColor = nitriteColor;
activeHColor = substrateColor;

labelTextSize = 14;
textSize = 12;
font = 'Arial';
marker = '-';
lineWidth = 2;

% preliminary plot of experimental data
figure(1);
set(1, 'Position', [152   492   637   329], 'PaperPositionMode', 'auto');
h = plot(timeExp, nh4Exp, 'kd'), hold on;
set(h, 'MarkerFaceColor', [0 0 0]);
h = plot(timeExp, no3Exp, 'k^');
set(h, 'MarkerFaceColor', [0 0 0]);
h = plot(timeExp, po4Exp, 'ko');
%h = plot(timeExp, acExp, 'kd');

%% get the results

results = getResultsFromDirectory(dirName);

%% plot the simulation cycle data
figure(1);

% create a command to format each of the subplots
plotFormat = ['set(gca,''Color'', [1 1 1]),''FontSize'', 8'];
plotAxis = ['set(gca,''XLim'', [0 730])'];

% recover the variables
iteration = results.iteration.variable(1).value;
time = results.iteration.variable(2).value;
realTime = results.iteration.variable(3).value;
biovolume = results.iteration.variable(4).value;
runLengthX = results.iteration.variable(5).value;
runLengthY = results.iteration.variable(6).value;
runLengthZ = results.iteration.variable(7).value;
detachedBiomass = results.iteration.variable(8).value;
erodedBiomass = results.iteration.variable(9).value;
sloughedBiomass = results.iteration.variable(10).value;
producedBiomass = results.iteration.variable(11).value;
biomass = results.iteration.variable(12).value;
oxygen = results.iteration.variable(13).value;
oxygenRate = results.iteration.variable(14).value;
ammonium = results.iteration.variable(15).value;
ammoniumRate = results.iteration.variable(16).value;
nitrite = results.iteration.variable(17).value;
nitriteRate = results.iteration.variable(18).value;
nitrate = results.iteration.variable(19).value;
nitrateRate = results.iteration.variable(20).value;
substrate = results.iteration.variable(21).value;
substrateRate = results.iteration.variable(22).value;
phosphate = results.iteration.variable(23).value;
phosphateRate = results.iteration.variable(24).value;
activeNH = results.iteration.variable(25).value;
inertNH = results.iteration.variable(26).value;
activeNO = results.iteration.variable(27).value;
inertNO = results.iteration.variable(28).value;
activeH = results.iteration.variable(29).value;
inertH = results.iteration.variable(30).value;
activePAO = results.iteration.variable(31).value;
phbPAO = results.iteration.variable(32).value;
polypPAO = results.iteration.variable(33).value;
glycogenPAO = results.iteration.variable(34).value;
inertPAO = results.iteration.variable(35).value;

totalInert = inertNH + inertNO + inertH + inertPAO;
timeDay = time/24;
averageRunLength = (runLengthX + runLengthY) * 0.5;


iterationToShowIndex = unique(find(mod(time, 15) == 0));

% find the indexes for the begining and end of the cycles
beginCycleTimes = 0:3:time(end);
[uniqueTime, iUT, jUT] = unique(time);

% cycle indexes:
%end indexes
estimatedEndIndexes =...
    interp1(uniqueTime, iUT, beginCycleTimes(2:end) - 0.001, 'linear') ;

endOfCycleIndex = round(estimatedEndIndexes);

%begining indexes
estimatedBeginIndexes =...
    interp1(uniqueTime, iUT, beginCycleTimes + 0.001, 'linear') ;

beginCycleIndex = ceil(estimatedBeginIndexes);
beginCycleIndex = beginCycleIndex(1:end-1);

%draw the plots
iterationIndex = length(iterationToShowIndex)-1;
ia = 1:iterationToShowIndex(iterationIndex); %index array
%find the time of the index of the end of the present
%cycle
interp1(uniqueTime, iUT, time(ia(end)) + 3, 'nearest');
ia2 = ia(end):interp1(uniqueTime, iUT, time(ia(end)) + 3, 'nearest');
% time in end of cycle array
ia3 =...
    endOfCycleIndex(find(time(endOfCycleIndex)<=time(ia(end))+3));

% plot compositions during cycle
matureStateIndexes = find(time > 150*24);

beginCycleIndexMature =...
    beginCycleIndex(beginCycleIndex > matureStateIndexes(1));
endOfCycleIndexMature =...
    endOfCycleIndex(endOfCycleIndex > matureStateIndexes(1));

% Compile data from mature cycles
timeMatureCycles = [];
substrateMatureCycles = [];
phosphateMatureCycles = [];
ammoniumMatureCycles = [];
nitriteMatureCycles = [];
nitrateMatureCycles = [];
oxygenMatureCycles = [];
phbPAOMatureCycles = [];
polypPAOMatureCycles = [];
glycogenPAOMatureCycles = [];
%for i = 400:400,
%for i = 1:100:(length(beginCycleIndexMature)-1),
for i = (length(beginCycleIndexMature)-1); 
    a = beginCycleIndexMature(i);
    b = endOfCycleIndexMature(i+1);
    c = a:b;
    timeCycle = [time(c) - time(a)];
    timeMatureCycles = [timeMatureCycles timeCycle'];
    substrateMatureCycles = [substrateMatureCycles substrate(c)'];
    phosphateMatureCycles = [phosphateMatureCycles phosphate(c)'];
    ammoniumMatureCycles = [ammoniumMatureCycles ammonium(c)'];
    nitriteMatureCycles = [nitriteMatureCycles nitrite(c)'];
    nitrateMatureCycles = [nitrateMatureCycles nitrate(c)'];
    oxygenMatureCycles = [oxygenMatureCycles oxygen(c)'];
    % and storage compounds
    phbPAOMatureCycles = [phbPAOMatureCycles phbPAO(c)'];
    polypPAOMatureCycles = [polypPAOMatureCycles polypPAO(c)'];
    glycogenPAOMatureCycles = [glycogenPAOMatureCycles glycogenPAO(c)'];
end;


anarobic = find (timeMatureCycles < 1);
aerobic = find (timeMatureCycles >= 1);
% itermittent lines
h = plot(timeMatureCycles(anarobic), phosphateMatureCycles(anarobic), marker);
set(h, 'Color', [0 0 0], 'LineWidth', lineWidth, 'LineStyle', '--');
hold on;
h = plot(timeMatureCycles(anarobic), nitrateMatureCycles(anarobic), marker);
set(h, 'Color', [0 0 0], 'LineWidth', lineWidth, 'LineStyle', '--');
h = plot(timeMatureCycles(anarobic), ammoniumMatureCycles(anarobic), marker);
set(h, 'Color', [0 0 0], 'LineWidth', lineWidth, 'LineStyle', '--');

% continuum lines
h = plot(timeMatureCycles(aerobic), phosphateMatureCycles(aerobic), marker);
set(h, 'Color', [0 0 0], 'LineWidth', lineWidth);
hold on;
h = plot(timeMatureCycles(aerobic), nitrateMatureCycles(aerobic), marker);
set(h, 'Color', [0 0 0], 'LineWidth', lineWidth);
h = plot(timeMatureCycles(aerobic), ammoniumMatureCycles(aerobic), marker);
set(h, 'Color', [0 0 0], 'LineWidth', lineWidth);
hold off;


xlabel('Time in cycle [h]', 'FontSize', labelTextSize, 'FontName', font);
ylabel({'Solutes in bulk' '[gN.L^{-1}, gP.L^{-1}]'},...
    'FontSize', labelTextSize, 'FontName', font);
set(gca, 'FontSize', textSize, 'FontName', font);
set(gca, 'FontSize', textSize, 'FontName', font);
set(gca, 'YLim', [0 0.12]);

print('-depsc', [resultPath '\temp\figure6-matlab.eps']);