%%% Enhances the annotation of a draft model on the basis of a similar,
%%% more accurate model and a list of ortholog genes between the 2

function enhancedModel = enhanceDraftModel(draftModel,topNotchModel,dictionary,gambleFlag)

%%% Added gambleFlag to control the third enhancement passage. If
% gambleFlag = 0, no excludeRxns will be put back into the model, as
% they've been proven to be a bit of a gamble (they're probably wrong
% anyway)

gambleFlag = logical(gambleFlag);

% consider "Unknown" and "Spontaneous" genes as universal
dictionary = cat(1,dictionary,{'Unknown','Unknown';'Spontaneous','Spontaneous'});

% Replace the _c0 and _e0 tag at the end of met names with a bracketed tag 
% (cobra doesn't like otherwise)
topNotchModel.mets = strrep(topNotchModel.mets,'_c0','[c]');
topNotchModel.mets = strrep(topNotchModel.mets,'_e0','[e]');
draftModel.mets = strrep(draftModel.mets,'_c0','[c]');
draftModel.mets = strrep(draftModel.mets,'_e0','[e]');

%Insert safety 'X's at the end of every gene, everywhere
for n=1:length(dictionary)
    dictionary{n,1}=regexprep(dictionary{n,1},'(?<gene>.*)','$<gene>X');
    dictionary{n,2}=regexprep(dictionary{n,2},'(?<gene>.*)','$<gene>X');
end
for n=1:length(topNotchModel.genes)
    topNotchModel.genes{n}=regexprep(topNotchModel.genes{n},'(?<gene>.*)','$<gene>X');
end
for n=1:length(topNotchModel.grRules)
    topNotchModel.grRules{n}=regexprep(topNotchModel.grRules{n},'(sm[^ )]*)','$1X');
end
for n=1:length(topNotchModel.grRules)
    topNotchModel.grRules{n}=strrep(topNotchModel.grRules{n},'Unknown[^X]','UnknownX');
end
for n=1:length(draftModel.genes)
    draftModel.genes{n}=regexprep(draftModel.genes{n},'(?<gene>.*)','$<gene>X');
end
for n=1:length(draftModel.grRules)
    draftModel.grRules{n}=regexprep(draftModel.grRules{n},'(An[^ )]*)','$1X');
end
for n=1:length(draftModel.grRules)
    draftModel.grRules{n}=strrep(draftModel.grRules{n},'Unknown[^X]','UnknownX');
end
for n=1:length(topNotchModel.grRules)
    topNotchModel.grRules{n}=strrep(topNotchModel.grRules{n},'Spontaneous[^X]','SpontaneousX');
end

% find out what common reactions exist between the draft and the topNotch
% model, and solve eventual name discrepancies
commonRxns = findCommonRxns(draftModel,topNotchModel);
for n = 1:length(commonRxns(:,1))
    if ~ismember(commonRxns(n,2),draftModel.rxns)
        tmpID = strmatch(commonRxns(n,1),draftModel.rxns,'exact');
        draftModel.rxns(tmpID) = strrep(draftModel.rxns(tmpID),commonRxns(n,1),commonRxns(n,2));
    end
end

%Add a .rxnEquations field to the 2 submodels
draftModel_rxnEquations = printRxnFormula(draftModel);
draftModel = setfield(draftModel,'rxnEquations',draftModel_rxnEquations);
topNotchModel_rxnEquations = printRxnFormula(topNotchModel);
topNotchModel = setfield(topNotchModel,'rxnEquations',topNotchModel_rxnEquations);

% identify the exchange and biomass reactions in the draft model
EXlist_draft = draftModel.rxns(strncmp('EX_',draftModel.rxns,3));
biomass_draft = draftModel.rxns(strncmp('biomass',draftModel.rxns,7));

% find out which reactions can be carried out by ortholog genes alone in 
% the accurate model and extract the relevant GPRs
[~,~,constrRxnNames,~] = deleteModelGenes(topNotchModel,setdiff(topNotchModel.genes,dictionary(:,1)));
sureRxns = setdiff(topNotchModel.rxns,constrRxnNames);
sureList = cell(length(sureRxns),2);
sureList(:,1) = sureRxns;
for n = 1:length(sureList(:,1))
    sureList(n,2) = topNotchModel.grRules(strmatch(sureList(n,1),topNotchModel.rxns,'exact'));
end

% find out which reactions are associated with ortholog genes, but also
% need other genes to work (ortholog AND non-ortholog)
%[~,list] = findRxnsFromGenes(topNotchModel,intersect(topNotchModel.genes,dictionary(:,1)),0,1);
%needyRxns = setdiff(unique(list(:,1)),sureRxns);
%needyList = cell(length(needyRxns),2);
%needyList(:,1) = needyRxns;
%for n = 1:length(needyList(:,1))
%    needyList(n,2) = topNotchModel.grRules(strmatch(needyList(n,1),topNotchModel.rxns,'exact'));
%end

% find out what novelty reactions brings this fierce young draft model,
% carried by his unique genes
%[~,list] = findRxnsFromGenes(draftModel,setdiff(draftModel.genes,dictionary(:,2)),0,1);
%newRxns = unique(list(:,1));
%newList = cell(length(newRxns),2);
%newList(:,1) = newRxns;
%for n = 1:length(newList(:,1))
%    newList(n,2) = draftModel.grRules(strmatch(newList(n,1),draftModel.rxns,'exact'));
%end


%%% start defining the reaction rules for the new, draft-no-more model


%% First, find out which of the sureRxns are not already in the draft model.
% For them, the GPRs shall be those of the topNotch, after sderenating the
% non-ortholog genes
sureNewRxns = setdiff(sureRxns,draftModel.rxns);
sureNewList = cell(length(sureNewRxns),8);
sureNewList(:,1) = sureNewRxns;
for n = 1:length(sureNewList(:,1))
    tmpID = strmatch(sureNewList(n,1),topNotchModel.rxns,'exact');
    sureNewList(n,2) = topNotchModel.rxnNames(tmpID);
    sureNewList(n,3) = topNotchModel.rxnEquations(tmpID);
    sureNewList{n,4} = topNotchModel.rev(tmpID);
    sureNewList{n,5} = topNotchModel.lb(tmpID);
    sureNewList{n,6} = topNotchModel.ub(tmpID);
    sureNewList(n,7) = topNotchModel.subSystems(tmpID);
    sureNewList(n,8) = topNotchModel.grRules(tmpID);
end
nonOrthologsTopNotch = setdiff(topNotchModel.genes,dictionary(:,1));
for n = 1:length(nonOrthologsTopNotch)
    sureNewList(:,8) = strrep(sureNewList(:,8),nonOrthologsTopNotch(n),'False');
end

% The more brackets the better
sureNewList(:,8) = regexprep(sureNewList(:,8),'(.*)','\($1\)');
%%% SDERENATE %%%

for n = 1:length(sureNewList(:,8))
    if ~isempty(sureNewList{n,8})
        tmpGPR = sprintf('''%s''',sureNewList{n,8});
        [~,tmpGPR] = system(sprintf('python booleanSderenator.py %s',tmpGPR));
        sureNewList(n,8) = cellstr(tmpGPR);
    end
end

% Too many brackets
sureNewList(:,8) = regexprep(sureNewList(:,8),'\((.*)\)','$1');

% Translate all in terms of draft orthologs
for n = 1:length(dictionary(:,1))
    sureNewList(:,8) = strrep(sureNewList(:,8),dictionary(n,1),dictionary(n,2));
end


%% Secondly, find out which of the sureRxns are already present in the draft
% model. Their fate is rather more complex. The final GPRs will be those of
% the topNotch, sderenated of any gene that has no ortholog, concatenated
% with an 'or' statement with the GPRs of the draft, the latter having been
% sderenated of any ortholog genes
sureOldRxns = intersect(sureRxns,draftModel.rxns);
sureOldList = cell(length(sureOldRxns),9);
sureOldList(:,1) = sureOldRxns;
for n = 1:length(sureOldList(:,1))
    tmpID = strmatch(sureOldList(n,1),topNotchModel.rxns,'exact');
    sureOldList(n,2) = topNotchModel.rxnNames(tmpID);
    sureOldList(n,3) = topNotchModel.rxnEquations(tmpID);
    sureOldList{n,4} = topNotchModel.rev(tmpID);
    sureOldList{n,5} = topNotchModel.lb(tmpID);
    sureOldList{n,6} = topNotchModel.ub(tmpID);
    sureOldList(n,7) = topNotchModel.subSystems(tmpID);
    sureOldList(n,8) = topNotchModel.grRules(tmpID);
    sureOldList(n,9) = draftModel.grRules(strmatch(sureOldList(n,1),draftModel.rxns,'exact'));
end
for n = 1:length(nonOrthologsTopNotch)
    sureOldList(:,8) = strrep(sureOldList(:,8),nonOrthologsTopNotch(n),'False');
end
orthologsDraft = intersect(draftModel.genes,dictionary(:,2));
for n = 1:length(orthologsDraft)
    sureOldList(:,9) = strrep(sureOldList(:,9),orthologsDraft(n),'False');
end

% The more brackets the better
sureOldList(:,8) = regexprep(sureOldList(:,8),'(.*)','\($1\)');
sureOldList(:,9) = regexprep(sureOldList(:,9),'(.*)','\($1\)');

%%% SDERENATE %%%
for n = 1:length(sureOldList(:,8))
    if ~isempty(sureOldList{n,8})
        tmpGPR = sprintf('''%s''',sureOldList{n,8});
        [~,tmpGPR] = system(sprintf('python booleanSderenator.py %s',tmpGPR));
        sureOldList(n,8) = cellstr(tmpGPR);
    end
end
for n = 1:length(sureOldList(:,9))
    if ~isempty(sureOldList{n,9})
        tmpGPR = sprintf('''%s''',sureOldList{n,9});
        [~,tmpGPR] = system(sprintf('python booleanSderenator.py %s',tmpGPR));
        sureOldList(n,9) = cellstr(tmpGPR);
    end
end

% Too many brackets
sureOldList(:,8) = regexprep(sureOldList(:,8),'\((.*)\)','$1');
sureOldList(:,9) = regexprep(sureOldList(:,9),'\((.*)\)','$1');

% Concatenate 8 and 9
for n = 1:length(sureOldList(:,1))
    if ~isempty(sureOldList{n,9})
        sureOldList(n,8) = strcat({'('},sureOldList(n,8),{' or '},sureOldList(n,9),{')'});
    end
end
sureOldList(:,9) = [];

%%% It may happen that, in the draft model, the only gene bringing the
%%% reaction about is an ortholog, and the first sderenation results in a
%%% 'False' statement. So a second sderenation is necessary.

% First, make sure there aren't any ' or 's floating around
sureOldList(:,8) = regexprep(sureOldList(:,8),'((?<=\() or | or (?=\)))','');

% The more brackets the better
sureOldList(:,8) = regexprep(sureOldList(:,8),'(.*)','\($1\)');

%%% SDERENATE %%%
for n = 1:length(sureOldList(:,8))
    if ~isempty(sureOldList{n,8})
        tmpGPR = sprintf('''%s''',sureOldList{n,8});
        [~,tmpGPR] = system(sprintf('python booleanSderenator.py %s',tmpGPR));
        sureOldList(n,8) = cellstr(tmpGPR);
    end
end

% Too many brackets
sureOldList(:,8) = regexprep(sureOldList(:,8),'\((.*)\)','$1');

% Translate all in terms of draft orthologs
for n = 1:length(dictionary(:,1))
    sureOldList(:,8) = strrep(sureOldList(:,8),dictionary(n,1),dictionary(n,2));
end

%% Thirdly (?), we have to make sure to exclude from the final model any
% reaction that is carried by ortholog genes alone, but is only present in
% the draft model, as we assume that if this reaction is not present in the
% topNotch model (while it could), then it shouldn't even be in the draft.
% Among these reactions, however, we should keep those which keep
% functioning thanks to some other gene when the ortholog is taken out of
% the relevant GPR. The latter (if there are any) should be kept after being contextually
% sderenated
[~,~,constrRxnNames] = deleteModelGenes(draftModel,setdiff(draftModel.genes,dictionary(:,2)));
ortOnlyRxns = setdiff(draftModel.rxns,constrRxnNames);
excludeRxns = setdiff(ortOnlyRxns,topNotchModel.rxns);
[~,~,constrRxnNames] = deleteModelGenes(draftModel,orthologsDraft);
nonOrtOnlyRxns = setdiff(draftModel.rxns,constrRxnNames);

if gambleFlag
    notExcludeRxns = setdiff(intersect(excludeRxns,nonOrtOnlyRxns),cat(1,biomass_draft,EXlist_draft));
else notExcludeRxns = [];
end

notExcludeList = cell(length(notExcludeRxns),3);

if ~isempty(notExcludeRxns)

notExcludeList(:,1) = notExcludeRxns;
for n = 1:length(notExcludeList(:,1))
    tmpID = strmatch(notExcludeList(n,1),draftModel.rxns,'exact');
    notExcludeList(n,2) = draftModel.rxnNames(tmpID);
    notExcludeList(n,3) = draftModel.rxnEquations(tmpID);
    notExcludeList{n,4} = draftModel.rev(tmpID);
    notExcludeList{n,5} = draftModel.lb(tmpID);
    notExcludeList{n,6} = draftModel.ub(tmpID);
    notExcludeList(n,7) = draftModel.subSystems(tmpID);
    notExcludeList(n,8) = draftModel.grRules(tmpID);
end
for n = 1:length(orthologsDraft)
    notExcludeList(:,8) = strrep(notExcludeList(:,8),orthologsDraft(n),'False');
end

% The more brackets the better (SAREBBE MEGLIO METTERLO DOPO IF, SIA QUA CHE ALTROVE)
notExcludeList(:,8) = regexprep(notExcludeList(:,8),'(.*)','\($1\)');
%%% SDERENATE %%%

for n = 1:length(notExcludeList(:,8))
    if ~isempty(notExcludeList{n,8})
        tmpGPR = sprintf('''%s''',notExcludeList{n,8});
        [~,tmpGPR] = system(sprintf('python booleanSderenator.py %s',tmpGPR));
        notExcludeList(n,8) = cellstr(tmpGPR);
    end
end

% Too many brackets
notExcludeList(:,8) = regexprep(notExcludeList(:,8),'\((.*)\)','$1');
end

% Finally, all the reactions in the draft model which do not belong in
% any of the previously defined categories (i.e. they are present in both
% models but with completely different GPRs, or they are exclusive to the 
% draft model) will be left as they are in the draft model, with their 
% current associations. But, in case the reaction is present in both of the
% models, upper and lower bounds and reversibility will be taken from the
% topNotch model.
% (A different case should be made for reactions bearing an 'and'
% relationship between a gene that is present in both models and something
% else that is not. But since there aren't any reactions common to both
% models in my current work that fall within this category, I'll just skip
% this step for now)
brandOldRxns = setdiff(draftModel.rxns,cat(1,sureOldRxns,excludeRxns,EXlist_draft,biomass_draft));

brandOldList = cell(length(brandOldRxns),3);
brandOldList(:,1) = brandOldRxns;
for n = 1:length(brandOldList(:,1))
    tmpID = strmatch(brandOldList(n,1),draftModel.rxns,'exact');
    brandOldList(n,2) = draftModel.rxnNames(tmpID);
    brandOldList(n,3) = draftModel.rxnEquations(tmpID);
    if ismember(brandOldList(n,1),topNotchModel.rxns)
        tmpID2 = strmatch(brandOldList(n,1),topNotchModel.rxns,'exact');
        brandOldList{n,4} = topNotchModel.rev(tmpID2);
        brandOldList{n,5} = topNotchModel.lb(tmpID2);
        brandOldList{n,6} = topNotchModel.ub(tmpID2);
    else brandOldList{n,4} = draftModel.rev(tmpID);
        brandOldList{n,5} = draftModel.lb(tmpID);
        brandOldList{n,6} = draftModel.ub(tmpID);
    end
    brandOldList(n,7) = draftModel.subSystems(tmpID);
    brandOldList(n,8) = draftModel.grRules(tmpID);
end

% Put all of the reactions in a single list and generate the freaking model
% already
if ~isempty(notExcludeList)
    enhancedList = cat(1,sureNewList,sureOldList,brandOldList,notExcludeList);
else enhancedList = cat(1,sureNewList,sureOldList,brandOldList);
end
rxnAbrList = enhancedList(:,1);
rxnNameList = enhancedList(:,2);
rxnList = enhancedList(:,3);
revFlagList = cell2mat(enhancedList(:,4));
lowerBoundList = cell2mat(enhancedList(:,5));
upperBoundList = cell2mat(enhancedList(:,6));
subSystemList = enhancedList(:,7);
grRuleList = enhancedList(:,8);

bigModel = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList,...
     lowerBoundList,upperBoundList,subSystemList,grRuleList,'','');

%Remove safety 'X's
bigModel.grRules = strrep(bigModel.grRules,'X','');
bigModel.genes = strrep(bigModel.genes,'X','');

enhancedModel = bigModel;

end