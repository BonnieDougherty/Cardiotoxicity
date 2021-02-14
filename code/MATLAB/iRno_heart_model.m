%% Generating rat iCardio model from iRno and iCardio
clear all
close all
clc

% Load the COBRA toolbox
initCobraToolbox
changeCobraSolver('ibm_cplex', 'all');

% Loading the iRno COBRA file
rno_cobra_load = ncomm_blais_xls2model('data/ncomms14250-s12, iRno COBRA.xlsx');
% Change objective function for ATP synthesis
rno_cobra_load = changeObjective(rno_cobra_load, 'RCR11017');

% Loading the iHsa COBRA file
hsa_cobra_load = ncomm_blais_xls2model('data/ncomms14250-s10, iHsa COBRA.xlsx');
% Change objective function for ATP synthesis
hsa_cobra_load = changeObjective(hsa_cobra_load, 'RCR11017');

% Make changes to the original reconstruction based off new metabolic tasks
% 1: Change the DHAP reaction to be reveresed
f = findRxnIDs(rno_cobra_load,'RCR21050');
rno_cobra_load.lb(f) = -1000;
rno_cobra_load.ub(f) = 0;
% Make the same changes for hsa_cobra_load
f = findRxnIDs(hsa_cobra_load,'RCR21050');
hsa_cobra_load.lb(f) = -1000;
hsa_cobra_load.ub(f) = 0;

% 2: Remove superoxide movement from the model (currently just exported out)
f = findRxnIDs(rno_cobra_load,'RCR40428');
rno_cobra_load.lb(f) = 0;
rno_cobra_load.ub(f) = 0;
% Make the same changes for hsa_cobra_load
f = findRxnIDs(hsa_cobra_load,'RCR40428');
hsa_cobra_load.lb(f) = 0;
hsa_cobra_load.ub(f) = 0;

% 3: All catalase reactions have bounds of [0,0] change to [0,1000]
rxns = findRxnIDs(rno_cobra_load, {'RCR10165' 'RCR10607' 'RCR11007' 'RCR11029' 'RCR14124'});
rno_cobra_load.lb(rxns) = 0;
rno_cobra_load.ub(rxns) = 1000;
% Make the same changes for hsa_cobra_load
hsa_cobra_load.lb(rxns) = 0;
hsa_cobra_load.ub(rxns) = 1000;

% 4: Change the Complex I reaction to include superoxide production
rno_cobra_load = addReaction(rno_cobra_load, 'RCR21048', {'m03103[m]','m02553[m]','m02039[m]', 'm02630[m]', 'm03102[m]','m02552[m]','m02039[c]','m02631[m]'}, [-1 -1 -5 -0.04 1 1 4 0.04], false);
hsa_cobra_load = addReaction(hsa_cobra_load, 'RCR21048', {'m03103[m]','m02553[m]','m02039[m]', 'm02630[m]', 'm03102[m]','m02552[m]','m02039[c]','m02631[m]'}, [-1 -1 -5 -0.04 1 1 4 0.04], false);

% 5: Change the number of protons moved to generate one ATP
rno_cobra_load = addReaction(rno_cobra_load, 'RCR20085', {'m02751[m]','m01285[m]','m02039[c]', 'm01371[m]', 'm02039[m]','m02040[m]'}, [-1 -1 -2.7 1 2.7 1], true);
hsa_cobra_load = addReaction(hsa_cobra_load, 'RCR20085', {'m02751[m]','m01285[m]','m02039[c]', 'm01371[m]', 'm02039[m]','m02040[m]'}, [-1 -1 -2.7 1 2.7 1], true);

% 6: Change all exchange reaction bounds that are not 1000, 0, or -1000
lbx = find(rno_cobra_load.lb ~= -1000 & rno_cobra_load.lb ~= 0);
rno_cobra_load.lb(lbx) = -1000;
ubx = find(rno_cobra_load.ub ~= 1000 & rno_cobra_load.ub ~= 0);
rno_cobra_load.ub(ubx) = 1000;
% Make the same changes for the iHsa model
lbx = find(hsa_cobra_load.lb ~= -1000 & hsa_cobra_load.lb ~= 0);
hsa_cobra_load.lb(lbx) = -1000;
ubx = find(hsa_cobra_load.ub ~= 1000 & hsa_cobra_load.ub ~= 0);
hsa_cobra_load.ub(ubx) = 1000;

%% Compare iRno to iCardio to develop rat-specific cardiac model
% Load in the heart model
% stored as the variable heart_model_curation
load('data/HeartModel.mat')

% deleted contains rxn_ids for reactions removed from the model
deleted = [];
for i = 1:length(hsa_cobra_load.rxns)
    temp = strcmp(hsa_cobra_load.rxns{i}, heart_model_curation.rxns());
    if max(temp) == 0
        deleted(end+1,1) = i;
    end
end

% Re-create heart model based off new deleted reaction list
deleted_rxns = {};
for k = 1:length(deleted)
    deleted_rxns{end+1} = hsa_cobra_load.rxns{deleted(k)};
end

rno_heart_model = removeRxns(rno_cobra_load, deleted_rxns);
hsa_heart_model = removeRxns(hsa_cobra_load, deleted_rxns);

% Remove gluconeogenesis from the model through transport of glucose out of
% the model
temp = findRxnIDs(rno_heart_model,'RCR40464');
rno_heart_model.ub(temp) = 0;

%% Check metabolic tasks with the new rat heart model
model = rno_heart_model;

% Check completion of metabolic tasks from iCardio
inputFile = ['data/AllTasks_CardiomyocyteSpecific_RAT.xlsx'];
[FINAL] = generateCobraTaskList(inputFile, rno_cobra_load);
xlswrite('data/AllTasks_CardiomyocyteSpecific_RAT_COBRA.xlsx', FINAL, 'TASKS')
inputFile = ['data/AllTasks_CardiomyocyteSpecific_RAT_COBRA.xlsx'];

rno_tasks = checkMetabolicTasks_BVD(model,inputFile);

%% Look at the reactions that are in iRno that aren't in iHsa
% All reactions are in both models, bounds are different
% 72 reactions are disabled in iHsa
iHsa_removed = [];
for rxn = 1:length(hsa_cobra_load.rxns)
    if(hsa_cobra_load.lb(rxn) == 0 && hsa_cobra_load.ub(rxn) == 0)
        iHsa_removed(end+1,1) = rxn;
    end
end

% Only 61 reactions disabled in iRno
iRno_removed = [];
for rxn = 1:length(rno_cobra_load.rxns)
    if(rno_cobra_load.lb(rxn) == 0 && rno_cobra_load.ub(rxn) == 0)
        iRno_removed(end+1,1) = rxn;
    end
end

iRno_removed = rno_cobra_load.rxns(iRno_removed);
iHsa_removed = hsa_cobra_load.rxns(iHsa_removed);

% Compare to see the 11 reactions that are different
for rxn = 1:length(iHsa_removed)
    if sum(strcmp(iHsa_removed{rxn}, iRno_removed)) == 0
        iHsa_removed{rxn}
    end
end

%% Reactions to be added based on NCBI evidence
add_rxns = {
    'RCR10020'
    'RCR90016'
    'RCR90127'
    'RCR90128'
    'RCR90129'
    'RCR90135'
    'RCR90145'
    'RCR90148'
    % no evidence of expression in the heart but included all reactions
    % associated with this GPR
    'RCR10559'
    'RCR11533'
    'RCR14431'
    'RCR90004'
    'RCR90005'};

% deleted contains rxn_ids for reactions removed from the model
deleted = [];
for i = 1:length(rno_cobra_load.rxns)
    temp = strcmp(rno_cobra_load.rxns{i}, rno_heart_model.rxns());
    if max(temp) == 0
        deleted(end+1,1) = i;
    end
end

deleted_update = deleted;
for i = 1:length(add_rxns)
    [trash rxn_ID] = max(strcmp(add_rxns{i}, rno_cobra_load.rxns()));
    location = find(deleted_update == rxn_ID);
    deleted_update(location) = [];
end

% Re-create heart model based off new deleted reaction list
deleted_rxns = {};
for k = 1:length(deleted_update)
    deleted_rxns{end+1} = rno_cobra_load.rxns{deleted_update(k)};
end

rno_heart_model = removeRxns(rno_cobra_load, deleted_rxns);

% Remove gluconeogenesis from the model through transport of glucose out of
% the model
temp = findRxnIDs(rno_heart_model,'RCR40464');
rno_heart_model.ub(temp) = 0;

%% Re-check metabolic tasks
model = rno_heart_model;
rno_tasks = checkMetabolicTasks_BVD(model,inputFile);
%% Save the model to a new MATLAB file
save('data/iRno_heart.mat', 'rno_heart_model', 'rno_cobra_load')