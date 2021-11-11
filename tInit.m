function tInit(tissue)
tissue
tpm_data = readtable("counts/TPMaverage.csv");
tpm_data(1:5, 1:5)


% extract the tissue and gene names
data_struct.tissues = tpm_data.Properties.VariableNames(2:end)';  % sample (tissue) names
data_struct.genes = tpm_data.Var1;  % gene names
data_struct.levels = table2array(tpm_data(:, 2:end));  % gene TPM values
data_struct.threshold = 1;

data_struct

load("models/Mouse-GEM.mat")


mouseGEM = addBoundaryMets(mouseGEM)

essentialTasks = parseTaskList('~/NoApps/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx')

checkTasks(mouseGEM, [], true, false, false, essentialTasks);

refModel = mouseGEM;
%tissue = 'ob_ob_HFD_Ao';
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = data_struct;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params = [];  % additional optimization parameters for the INIT algorithm
paramsFT = [];  % additional optimization parameters for the fit-tasks algorithm
params.TimeLimit = 5000;  % Increase timelimits


model= getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

ModelName = strcat('model_', tissue);
eval([strcat(ModelName, '= model')])

% export sbml models
exportModel(eval(ModelName),strcat('models/context/', ModelName,'.xml'))
% export matlab models
save(strcat('models/context/', ModelName,'.mat'), ModelName)
end
