function results = getResultsFromDirectory(dirName)

% results = getResultsFromDirectory(dirName)
% Parse directory structure for results

results.directory = dirName;
results.validData = false;
results.particles = [];
results.iteration = [];
results.solutes   = [];
results.solids    = [];

% check if all results subdirectories are available
if ~(isdir([dirName '/particles']) ...
        & isdir([dirName '/solids']) ...
        & isdir([dirName '/solutes']))
    return;
end;

% get the iteration numbers from the particle file listings
d = dir([dirName '/particles/iteration*']);
preNumber = length('iteration');
postNumber = length('.txt');
% return if no data is found
if (length(d) == 0)
    return;
end;

% get the iteration numbers and the times from simulationParameters.txt
results = readSimulationParameters([dirName '/simulationResults.txt'], results);
results.iteration.current = 2;

% set default value for the particle checkbox
results.particles = true;

%SOLUTES
% setup up default value of solute to true (show solute concentration)
results.solutes.show = true;
% get the names of solutes from the first iteration
d = dir([dirName '/solutes/iteration00000000/*.txt']);
for i = 1:length(d),
    results.solutes.available{i} = d(i).name(1:end-postNumber);
end;
results.solutes.current = 1;

%SOLIDS
% setup up default value of solids to false (show solute concentration)
results.solids.show = false;
% get the names of solutes from the first iteration
d = dir([dirName '/solids/iteration00000000/*.txt']);
for i = 1:length(d),
    results.solids.available{i} = d(i).name(1:end-postNumber);
end;
results.solids.current = 1;

%FLOW
% setup up default value of flow to false (show solute concentration)
results.flow.show = false;
% get the names of solutes from the first iteration
d = dir([dirName '/flow/iteration00000000/*.txt']);
for i = 1:length(d),
    results.flow.available{i} = d(i).name(1:end-postNumber);
end;
results.flow.current = 1;

% validate the structure
results.validData = true;
