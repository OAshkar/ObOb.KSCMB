DIR = "./models/context/";
contextfiles = dir(DIR);
contextfiles = regexpi({contextfiles.name}, ".*mat", "match", "once");
contextfiles = contextfiles(~cellfun('isempty',contextfiles)) % remove empty
contextfiles2 = regexpi(contextfiles,'model_.*(?=\.mat)','match')


cleannames = {} % only the names to be displayed
for i = 1:numel(contextfiles2)
    cleannames{i} = strrep(char(contextfiles2{i}), "_", ".");
    cleannames{i} = strrep(cleannames{i}, "model", "")
end


models = {}
for i = 1:numel(contextfiles)
    load(strcat(DIR, contextfiles{i}));
    %quick cleaning for IDs
    eval(strcat(char(contextfiles2{i}), " = setfield(", char(contextfiles2{i}), ",'id',", "'",cleannames{i}, "')"));
    models{i} = eval(char(contextfiles2{i}));
end

model_ids = arrayfun(@(i) models{i}.id, (1:numel(models))', 'UniformOutput', false);

res = compareMultipleModels(models); % compare the models
res

clustergram(res.structComp, 'Symmetric', false, 'Colormap', 'bone', 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);

writematrix(res.structComp,'MatResults/structComp.csv')
writecell(res.modelIDs,'MatResults/structComp_rowcolnames.csv')


rxn2Dmap = tsne(res.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 5);
% plot and label the GEMs in tSNE space
scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
hold on
text(rxn2Dmap(:,1), rxn2Dmap(:,2), res.modelIDs);
% TODO

subMat = res.subsystems.matrix;
subCoverage = (subMat - mean(subMat, 2)) ./ mean(subMat, 2) * 100;
% select subsystems to include in plot
inclSub = any(abs(subCoverage) > 25, 2);
subNames = res.subsystems.ID(inclSub);
% generate clustergram
cg = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap, 'DisplayRange', 100, 'rowLabels', subNames, 'columnLabels', [cleannames{:}], 'ShowDendrogram', 'OFF');

writematrix(subCoverage(inclSub,:),'MatResults/subCoverage.csv')
writecell(subNames,'MatResults/subCoverage_rownames.csv')
writematrix([cleannames{:}],'MatResults/subCoverage_colnames.csv')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare GEM functions

% Filter models
taskFile = '~/NoApps/Human-GEM/data/metabolicTasks/metabolicTasks_Full.xlsx' % from Human-GEM

res_func = compareMultipleModels(models, false, false, [], true, taskFile);
res_func.funcComp;

isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2);
diffTasks = res_func.funcComp.tasks(isDiff)

% visualize the matrix
spy(res_func.funcComp.matrix(isDiff,:), 30);

% apply some formatting changes
set(gca, 'XTick', 1:numel(cleannames), 'XTickLabel', [cleannames{:}], 'XTickLabelRotation', 90, 'YTick', 1:numel(diffTasks), 'YTickLabel', diffTasks, 'YAxisLocation', 'right');
xlabel(gca, '');

writematrix(res_func.funcComp.matrix(isDiff,:),'MatResults/MetabolicFunction.csv')
writecell(diffTasks,'MatResults/MetabolicFunction_rownames.csv')
writematrix([cleannames{:}],'MatResults/MetabolicFunction_colnames.csv')
