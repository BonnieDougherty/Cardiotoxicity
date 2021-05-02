% Testing various data sets with TIDEs
clear all
clc

% initialize the toolbox
initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

% Load in the heart model
load('data/iRno_heart.mat')

model = rno_heart_model;

% Generate a taskStructure for determining the min set of reactions
% Convert from original format to COBRA format
inputFile = ['data/AllTasks_CardiomyocyteSpecific_RAT_COBRA.xlsx'];
taskStructure=generateTaskStructure_BVD(inputFile);

% calculate the reactions necessary for each task
removeNoGPR = 'true';
minRxnList = generateMinRxnList(model, taskStructure, removeNoGPR);

% Add on analysis by reaction subsystem
% Find a unique list of subsystems in this model
reactions = {};
for k = 1:length(rno_heart_model.subSystems)
    reactions{k,1} = char(model.subSystems{k});
end
subsystem = unique(reactions);
subsystem(1,:) = [];

for k = 1:length(subsystem)
    minRxnList(end+1).id = strcat('S',num2str(k));
    minRxnList(end).description = subsystem{k};
    minRxnList(end).rxns = model.rxns(strcmp(subsystem{k}, reactions));
end

% Remove reactions from the subsystem list that don't have GPR rules
for k = 1:length(model.rxns)
    noGPR(k) = isempty(model.grRules{k});
end
noGPR = model.rxns(noGPR);

% Remove these reactions from the list of total_carrying_rxns
for task = 1:length(minRxnList)
    common = intersect(noGPR, minRxnList(task).rxns);
    minRxnList(task).rxns = setdiff(minRxnList(task).rxns, common);
end

% Remove tasks that no longer have associated reactions
task = 1;
while task ~= length(minRxnList)
    if length(minRxnList(task).rxns) < 3
        minRxnList(task) = [];
    else
        task = task+1;
    end
end

% Specify the number of random iterations of data
parpool(4);
num_iterations = 1000;

%% Debugging reading in table problems
load('opts.mat')
opts.Delimiter = {','};
opts.VariableNames = {'EntrezID','logfc','pval','fdr','ave'};

%% Run TIDEs analysis for cardiotoxicity data
% Load in data that had data for all genes, not just genes that are DEGs
% Read in data from RNA-seq analysis
Dox6hrs = readtable('data/RNA-seq/dougherty_rno_cardio_t6_dox_gene_deseq2.csv', opts);
Dox24hrs = readtable('data/RNA-seq/dougherty_rno_cardio_t24_dox_gene_deseq2.csv', opts);

Ace6hrs = readtable('data/RNA-seq/dougherty_rno_cardio_t6_ace_gene_deseq2.csv', opts);
Ace24hrs = readtable('data/RNA-seq/dougherty_rno_cardio_t24_ace_gene_deseq2.csv', opts);

FiveFU6hrs = readtable('data/RNA-seq/dougherty_rno_cardio_t6_5FU_gene_deseq2.csv', opts);
FiveFU24hrs = readtable('data/RNA-seq/dougherty_rno_cardio_t24_5FU_gene_deseq2.csv', opts);

%% Run test expression mapping to make sure EntrezIDs are formatted correctly
data.gene = Dox6hrs.EntrezID;
data.value = Dox6hrs.logfc;
[expressionRxns parsedGPR] = mapExpressionToReactions(model, data);
%%
% Run Dox6 hrs analysis
data.gene = Dox6hrs.EntrezID;
data.value = Dox6hrs.logfc;
[Dox6hrs, Dox6hrs_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Run Dox24 hrs analysis
data.gene = Dox24hrs.EntrezID;
data.value = Dox24hrs.logfc;
[Dox24hrs, Dox24hrs_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Run Ace6 hrs analysis
data.gene = Ace6hrs.EntrezID;
data.value = Ace6hrs.logfc;
[Ace6hrs, Ace6hrs_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Run Ace24 hrs analysis
data.gene = Ace24hrs.EntrezID;
data.value = Ace24hrs.logfc;
[Ace24hrs, Ace24hrs_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Run Five5FU 6 hrs analysis
data.gene = FiveFU6hrs.EntrezID;
data.value = FiveFU6hrs.logfc;
[FiveFU6hrs, FiveFU6hrs_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

% Run Five5FU 24 hrs analysis
data.gene = FiveFU24hrs.EntrezID;
data.value = FiveFU24hrs.logfc;
[Five5FU24hrs, Five5FU24hrs_random] = calculateTIDEscores(model, minRxnList, data, num_iterations);

%% Save all the datasets to Excel spreadsheets
data = Dox6hrs;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/TIDEs/Dox6hrs.xlsx', data_save)

data = Dox24hrs;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/TIDEs/Dox24hrs.xlsx', data_save)

%% Ace data
data = Ace6hrs;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/TIDEs/Ace6hrs.xlsx', data_save)

data = Ace24hrs;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/TIDEs/Ace24hrs.xlsx', data_save)

%% FiveFU
data = FiveFU6hrs;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/TIDEs/FiveFU6hrs.xlsx', data_save)

data = Five5FU24hrs;
% Save individual variable names to xlsx files in R/data folder
data_save = {};
for k = 1:length(data)
    data_save{k,2} = data(k).description;
    data_save{k,1} = data(k).id;
    data_save{k,4} = data(k).significance;
    data_save{k,3} = data(k).taskScore;
end
xlswrite('data/TIDEs/FiveFU24hrs.xlsx', data_save)

%% Save the random data sets
data = Ace6hrs_random; 
xlswrite('data/TIDEs/Ace6hrs_random.xlsx', data);

data = Dox24hrs_random;
xlswrite('data/TIDEs/Dox24hrs_random.xlsx', data);

data = Ace24hrs_random; 
xlswrite('data/TIDEs/Ace24hrs_random.xlsx', data);

data = Dox6hrs_random;
xlswrite('data/TIDEs/Dox6hrs_random.xlsx', data);

data = FiveFU6hrs_random; 
xlswrite('data/TIDEs/5FU6hrs_random.xlsx', data);

data = Five5FU24hrs_random;
xlswrite('data/TIDEs/5FU24hrs_random.xlsx', data);

%% Save the results to a MATLAB file
save('data/TIDEs/cardiotoxicity.mat','minRxnList','Ace24hrs','Ace6hrs','Dox24hrs','Dox6hrs','FiveFU24hrs','FiveFU6hrs')

%% Save list of reactions from metabolic tasks
ID = [];
reactions = [];
for task = 1:length(minRxnList)
    ID = [ID;repmat(string(minRxnList(task).id), length(minRxnList(task).rxns),1)];
    reactions = [reactions; minRxnList(task).rxns];
end

temp = table(ID, reactions);
writetable(temp, 'data/TIDEs/TIDEs_rxns.csv', 'Delimiter', ',')
