function Y = plotMappingInfo(resPath, patOrg, reacPat, reacTab, reacNumber, infoFilePath, figForm, sampNames, organisms)
% This function computes and automatically plots information coming from
% the mapping data as metabolic diversity and classical multidimensional
% scaling of individuals' reactions repertoire. If the last 2 arguments are 
% specified MDS plots will be annotated with samples and organisms names
%
% USAGE:
%
%   Y = plotMappingInfo(resPath, patOrg, reacPat, reacTab, reacNumber, infoFilePath, figForm, sampNames, organisms)
%
% INPUTS:
%   resPath:            char with path of directory where results are saved
%   reac:               nx1 cell array with all the unique set of reactions
%                       contained in the models
%   micRea:             binary matrix assessing presence of set of unique
%                       reactions for each of the microbes
%   reacSet:            matrix with names of reactions of each individual
%   reacTab:            binary matrix with presence/absence of reaction per
%                       individual.
%   reacAbun:           matrix with abundance of reaction per individual
%   reacNumber:         number of unique reactions of each individual
%   infoFilePath:       char indicating, if stratification criteria are available, 
%                       full path and name to related documentation(default: no)
%                       is available
%   figForm:            format to use for saving figures
%   sampNames:           nx1 cell array cell array with names of individuals in the study
%   organisms:          nx1 cell array cell array with names of organisms in the study
%
% OUTPUTS:
%   Y:                  classical multidimensional scaling of individuals'
%                       reactions repertoire
%
% .. Author: - Federico Baldini, 2017-2018

if ~exist('sampNames', 'var')
    sampNames = 0;
    aN = 0;
else
    aN=1;
end

if ~exist('organisms', 'var')
    organisms = 0;
    aO = 0;
else
    aO=1;
end

figure
imagesc(cell2mat(table2cell(reacPat(:,2:5))));
xlabel('Individuals');  % x-axis label
ylabel('Organisms');  % y-axis label
ax = gca;
ax.XTick = [1:length(patOrg)];
if  aO>0
    ax.YTick = [1:length(organisms)];
    organisms2=strrep(organisms,'_',' ');
    ax.YTickLabel = organisms2;
end
if  aN>0 
    ax.XTickLabel = sampNames;
    ax.XTickLabelRotation = 45;
end
title('Heatmap individuals | organisms reactions');
c = colorbar;
c.Label.String = 'Number of reactions';
print(strcat(resPath, 'Heatmap'), figForm)

if ~exist('infoFilePath', 'var')||~exist(infoFilePath, 'file')
    patStat = 0;
else
    patStat = 1;
end

if patStat == 0
% Plot:metabolic diversity
figure
A=[patOrg' reacNumber'];
[Auniq,~,IC] = unique(A,'rows');
cnt = accumarray(IC,1);

scatter(Auniq(:,1), Auniq(:,2), (10+5*(cnt>1)).^2); % make the ones where we'll put a number inside a bit bigger
for ii=1:numel(cnt)
    if cnt(ii)>1
        text(Auniq(ii,1),Auniq(ii,2),num2str(cnt(ii)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize', 6);
    end
end
ax = gca;
ax.XTick = [min(patOrg):max(patOrg)];
%xlim([min(patOrg) max(patOrg)]);
if length(unique(reacNumber))>1
    ylim([min(reacNumber) max(reacNumber)]);
end
xlabel('Microbiota Size')  % x-axis label
ylabel('Number of unique reactions')  % y-axis label
title('Metabolic Diversity')
print(strcat(resPath, 'Metabolic_Diversity'), figForm)

% PCoA -> different reactions per individual
D = pdist(cell2mat(table2cell(reacTab(:,2:end))),'jaccard');
[Y, eigvals] = cmdscale(D);
% catch if too few individuals to cluster
if ~isempty(Y)
    if (length(Y(1,:))>1)
        figure
        P = [eigvals eigvals / max(abs(eigvals))];
        expr = [eigvals/sum(eigvals)];
        plot(Y(:, 1), Y(:, 2), 'bx')
        title('PCoA of reaction presence');
        xlabel(strcat('PCoA1: ',num2str(round(expr(1)*100,2)),'% of explained variance'));
        ylabel(strcat('PCoA2: ',num2str(round(expr(2)*100,2)),'% of explained variance'));
%         if aN>0
%             text(Y(:,1),Y(:,2),sampNames,'HorizontalAlignment','left');%to insert numbers
%         end
        print(strcat(resPath, 'PCoA reactions'), figForm)
    else
        disp('noPcoA will be plotted')     
    end
end  

% build numbers of patients
% lab = 1:length(Y(:,1));
% lab = strread(num2str(a),'%s');
% labels = lab';
% text(Y(:,1),Y(:,2),labels,'HorizontalAlignment','left');%to insert numbers


else
    % Plot: number of species | number of reactions  disease resolved
    % Patients status: cellarray of same lenght of number of patients 0 means patient with disease 1 means helthy
infoFile = readInputTableForPipeline(infoFilePath);

% remove individuals not in simulations
[C,IA] = setdiff(infoFile(:,1),sampNames);
infoFile(IA,:)=[];

% get the number of conditions
cond=unique(infoFile(:,2));

% assign a random color to each condition

for i=1:length(cond)
    cols(i,:)=[rand rand rand];
end

colorMap = zeros(size(infoFile,1),3);

for i = 1: size(infoFile(:,1))
    % get the color corresponding to the condition
    findCol=find(strcmp(cond,infoFile{i,2}));
    colorMap(i, :) = cols(findCol,:);
end

figure
scatter(patOrg, reacNumber, 30 * ones(length(reacNumber), 1), colorMap, 'filled');
xlabel('Microbiota Size')  % x-axis label
ylabel('Number of unique reactions')  % y-axis label
title('Metabolic diversity')

for i=1:length(cond)
    text(max(patOrg(:, 1)-0.4*i),max(reacNumber(:, 2))+50,cond{i},'HorizontalAlignment','left','Color', cols(i,:));
        hold on
end

print(strcat(resPath, 'Metabolic_Diversity'), figForm)

% PCoA -> different reactions per individual
D = pdist(cell2mat(table2cell(reacTab(:,2:end))),'jaccard');
[Y, eigvals] = cmdscale(D);

P = [eigvals eigvals / max(abs(eigvals))];
if (length(Y(1,:))>2)
    figure
    expr = [eigvals/sum(eigvals)];
    plot(Y(:, 1), Y(:, 2), 'bx')
    title('PCoA of reaction presence');
    xlabel(strcat('PCoA1: ',num2str(round(expr(1)*100,2)),'% of explained variance'));
    ylabel(strcat('PCoA2: ',num2str(round(expr(2)*100,2)),'% of explained variance'));
    
%     for i=1:length(cond)
%     text(max(Y(:, 1)-0.02*i),max(Y(:, 2)),cond{i},'HorizontalAlignment','left','Color', cols(i,:));
%         hold on
%     end

    print(strcat(resPath, 'PCoA_reactions'), figForm)
else
    disp('noPcoA will be plotted')
end

end

% Plot Eigen number value: diasbled by default
% plot(1:length(eigvals),eigvals,'bo-');
% line([1,length(eigvals)],[0 0],'LineStyle',':','XLimInclude','off',...
%  'Color',[.7 .7 .7])
% axis([1,length(eigvals),min(eigvals),max(eigvals)*1.1]);
% xlabel('Eigenvalue number');
% ylabel('Eigenvalue');
% print(strcat(resPath,'Eigen number value'),figform)

% 3D PCoA plot
% scatter3(Y(:,1),Y(:,2),Y(:,3))
% print(strcat(resPath,'3D PCoA reactions'),figForm)
% text(Y(:,1),Y(:,2),Y(:,3),labels,'HorizontalAlignment','left');%to insert numbers

end
