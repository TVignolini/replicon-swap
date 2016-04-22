function [swappedModel,doubleRxns,duplicateRxns,redundantGenes] = repliconSwap(receivingModel,donorModel,recRemoveGenes,addRepliconGenes,dictionary) 

%Given a donor model, a set of replicon-associated genes for each replicon
%in the donor model, and a set of donor model genes which have at least an
%ortholog harbored by a receiving chromosome, this script generates a
%vector containing all the reactions to be transplanted in the receiving
%model along with a selected replicon or set of replicons

%Restrict the dictionary to only encompass genes already present in the
%receiving model
for n = 1:length(dictionary(:,1))
    if ~ismember(dictionary(n,1),receivingModel.genes)
        dictionary(n,:)={''};
    end
end
dictionary=dictionary(~cellfun('isempty',dictionary));
dictionary=reshape(dictionary',[length(dictionary)/2,2]);

%Replace the _c0 and _e0 tag at the end of met names with a bracketed tag (cobra
%doesn't like otherwise)
receivingModel.mets = strrep(receivingModel.mets,'_c0','[c]');
receivingModel.mets = strrep(receivingModel.mets,'_e0','[e]');
%donorModel.mets = strrep(donorModel.mets,'_c0','[c]');
%donorModel.mets = strrep(donorModel.mets,'_e0','[e]');

%Insert safety 'X's at the end of every gene, everywhere
%for n=1:length(badGenes)
%    badGenes{n}=regexprep(badGenes{n},'(?<gene>.*)','$<gene>X');
%end
for n=1:length(addRepliconGenes)
    addRepliconGenes{n}=regexprep(addRepliconGenes{n},'(?<gene>.*)','$<gene>X');
end
for n=1:length(recRemoveGenes)
    recRemoveGenes{n}=regexprep(recRemoveGenes{n},'(?<gene>.*)','$<gene>X');
end
for n=1:length(receivingModel.genes)
    receivingModel.genes{n}=regexprep(receivingModel.genes{n},'(?<gene>.*)','$<gene>X');
end
for n=1:length(receivingModel.grRules)
    receivingModel.grRules{n}=regexprep(receivingModel.grRules{n},'(sm[^ )]*)','$1X');
end
for n=1:length(receivingModel.grRules)
    receivingModel.grRules{n}=strrep(receivingModel.grRules{n},'Unknown','UnknownX');
end
for n=1:length(donorModel.genes)
    donorModel.genes{n}=regexprep(donorModel.genes{n},'(?<gene>.*)','$<gene>X');
end
for n=1:length(donorModel.grRules)
    donorModel.grRules{n}=regexprep(donorModel.grRules{n},'(An[^ )]*)','$1X');
end
for n=1:length(donorModel.grRules)
    donorModel.grRules{n}=strrep(donorModel.grRules{n},'Unknown','UnknownX');
end
for n=1:length(dictionary)
    dictionary{n,1}=regexprep(dictionary{n,1},'(?<gene>.*)','$<gene>X');
    dictionary{n,2}=regexprep(dictionary{n,2},'(?<gene>.*)','$<gene>X');
end

recModel=receivingModel;
donModel=donorModel;

%Create separate dictionaries
%dictionary_don2rec: dictionary used to translate non-transplanted genes in
%the donor model with receiving, chromosome-harbored orthologs
dictionary_don2rec=dictionary;
for n=1:length(dictionary)
    if ~isempty(strmatch(dictionary{n,2},addRepliconGenes,'exact'))
        dictionary_don2rec(n,:)={''};
    end
    if ~strncmp('smc',dictionary{n,1},3)
        dictionary_don2rec(n,:)={''};
    end
end
dictionary_don2rec=dictionary_don2rec(~cellfun('isempty',dictionary_don2rec));
dictionary_don2rec=reshape(dictionary_don2rec',[length(dictionary_don2rec)/2,2]);

%dictionary_rec2don: dictionary used to translate chromidic genes in the
%receiving model with incoming orthologs
dictionary_rec2don=dictionary;
for n=1:length(dictionary)
    if isempty(strmatch(dictionary{n,2},addRepliconGenes,'exact'))
        dictionary_rec2don(n,:)={''};
    end
    if strncmp(dictionary{n,1},'smc',3)
        dictionary_rec2don(n,:)={''};
    end
end
dictionary_rec2don=dictionary_rec2don(~cellfun('isempty',dictionary_rec2don));
dictionary_rec2don=reshape(dictionary_rec2don',[length(dictionary_rec2don)/2,2]);

%%% RECEIVING MODEL PREPARATION %%%

%Rename every gene to be deleted from the receiving model with the
%corresponding ortholog harbored by the incoming replicon or set of
%replicons
for n=1:length(dictionary_rec2don)
    disp((n/length(dictionary_rec2don))*100)
    recModel.grRules = strrep(recModel.grRules,dictionary_rec2don{n,1},dictionary_rec2don{n,2});
    recModel.genes = strrep(recModel.genes,dictionary_rec2don{n,1},dictionary_rec2don{n,2});
end

%If a gene in either of the models will end up in the swapped model along 
%with an ortholog, rename it as an 'OR' couple of the two orthologs in the 
%relevant GPRs
%First create a dictionary
dictionary_together=dictionary;
%Only consider genes that are going to be transplanted
for n=1:length(dictionary)
    if isempty(strmatch(dictionary{n,2},addRepliconGenes,'exact'))
        dictionary_together(n,:)={''};
    end
%Only consider chromosome-borne genes for the recipient
    if ~isempty(strmatch(dictionary{n,1},recRemoveGenes,'exact'))
        dictionary_together(n,:)={''};
    end
%Only consider genes that are at least in the receiving or in the donor model
    if (isempty(strmatch(dictionary{n,2},recModel.genes,'exact')) && isempty(strmatch(dictionary{n,2},donModel.genes,'exact')))
        dictionary_together(n,:)={''};
    end
end
dictionary_together=dictionary_together(~cellfun('isempty',dictionary_together));
dictionary_together=reshape(dictionary_together',[length(dictionary_together)/2,2]);
redundantGenes=strrep(dictionary_together,'X','');

%Generate a new set of grRules for both models
for n=1:length(dictionary_together(:,1))
    disp((n/length(dictionary_together(:,1)))*100)
    if ~isempty(strmatch(dictionary_together{n,2},donModel.genes,'exact'))
        geneID=strmatch(dictionary_together{n,2},donModel.genes,'exact');
        [~,list1]=findRxnsFromGenes(donModel,donModel.genes(geneID),'',1);
        for a=1:length(list1(:,1))
            rxnID1=strmatch(list1{a,1},donModel.rxns,'exact');
            newDonRule=strrep(donModel.grRules{rxnID1},dictionary_together{n,2},strcat({'('},dictionary_together(n,1),{' or '},dictionary_together(n,2),{')'}));
            donModel=changeGeneAssociation(donModel,donModel.rxns(rxnID1),cell2mat(newDonRule));
        end
        
    end
    if ~isempty(strmatch(dictionary_together{n,1},recModel.genes,'exact'))
        geneID=strmatch(dictionary_together{n,1},recModel.genes,'exact');
        [~,list2]=findRxnsFromGenes(recModel,recModel.genes(geneID),'',1);
        for b=1:length(list2(:,1))
            rxnID2=strmatch(list2{b,1},recModel.rxns,'exact');
            newRecRule=strrep(recModel.grRules{rxnID2},dictionary_together{n,1},strcat(dictionary_together(n,1),{' or '},dictionary_together(n,2)));
            recModel=changeGeneAssociation(recModel,recModel.rxns(rxnID2),cell2mat(newRecRule));
        end
    end
end        

%Removes all the reactions which are lost upon the removal of the chromids
%and creates a receiving submodel
[~, ~, constrRxnNames, ~] = deleteModelGenes(recModel,intersect(recModel.genes,recRemoveGenes));
Rmodel = extractSubNetwork(recModel,setdiff(recModel.rxns,constrRxnNames));

%%% DONOR MODEL PREPARATION %%%

%Remove unwanted reactions based on the curation of the original 1021 model
%if ortFlag
%    donModel = removeRxns(donModel,badReactions);
%end

%Rename every gene to be deleted from the donor model with the
%corresponding ortholog harbored by the recipient chromosome
for n=1:length(dictionary_don2rec)
    disp((n/length(dictionary_don2rec))*100)
    donModel.grRules = strrep(donModel.grRules,dictionary_don2rec{n,2},dictionary_don2rec{n,1});
    donModel.genes = strrep(donModel.genes,dictionary_don2rec{n,2},dictionary_don2rec{n,1});
end

%Simulate deletion of every non-transplanted gene in the donor model
%(except those having an ortholog in the receiving model) and select the
%reactions that keep working
[~, ~, constrRxnNames, ~] = deleteModelGenes(donModel,intersect(donModel.genes,setdiff(donorModel.genes,addRepliconGenes)));
[transplantRxns] = setdiff(donModel.rxns,constrRxnNames);

%Delete EX_ and 'biomass' reactions
transplantRxns = transplantRxns(strncmp('rxn',transplantRxns,3));

%Create raw submodel based on the reactions to be exported
Dmodel = extractSubNetwork(donModel,transplantRxns);

%%% PREVENTIVE SUBMODEL SDERENATION %%%

% Sderenate the submodels prior to fusion, so that the final swapped model
% comes out all tidy and clean

% REMOVE THE BAD GENES FROM THE RULES
unwantedGenesR = setdiff(intersect(Rmodel.genes,cat(1,recRemoveGenes,setdiff(donorModel.genes,addRepliconGenes))),'UnknownX');
%if ortFlag
%    badGenesR = intersect(Rmodel.genes,badGenes);
%    unwantedGenesR = unique(cat(1,unwantedGenesR,badGenesR));
%end
unwantedGenesD = setdiff(intersect(Dmodel.genes,cat(1,recRemoveGenes,setdiff(donorModel.genes,addRepliconGenes))),'UnknownX');
%if ortFlag
%    badGenesD = intersect(Dmodel.genes,badGenes);
%    unwantedGenesD = unique(cat(1,unwantedGenesD,badGenesD));
%end

% Rename the bad genes as 'False'
for n=1:length(unwantedGenesR)
    Rmodel.grRules = strrep(Rmodel.grRules,unwantedGenesR(n),'False');
end
for n=1:length(unwantedGenesD)
    Dmodel.grRules = strrep(Dmodel.grRules,unwantedGenesD(n),'False');
end

% Sderenate the submodels

Rmodel = sderenateModel(Rmodel);
Dmodel = sderenateModel(Dmodel);

%%% /PREVENTIVE SUBMODEL SDERENATION %%%

%Add the .rxnEquations fields to the 2 submodels
Dmodel_rxnEquations = printRxnFormula(Dmodel);
Dmodel = setfield(Dmodel,'rxnEquations',Dmodel_rxnEquations);
Rmodel_rxnEquations = printRxnFormula(Rmodel);
Rmodel = setfield(Rmodel,'rxnEquations',Rmodel_rxnEquations);

%Create a new model containing the sum of the submodels' genes and
%reactions and, in case of duplicate reactions, a combination of the 2 GPRs

%%% EDIT, APRIL 12th:
%%% The draft models are now devoid of any reaction abbreviation and gene
%%% association conflicts.
%%% In case of duplicate reactions, the gprs will be those of the receiving,
%%% connected with an ' or ' statement to those of the donor, the latter
%%% having been sderenated of any ortholog. It is likely that the final
%%% GPRs will only consist of the receiving ones.

doubleRxns = findCommonRxns(Rmodel,Dmodel);
wrongNames = setdiff(doubleRxns(:,2),doubleRxns(:,1));
rightRxns = unique(setdiff(cat(1,Rmodel.rxns,Dmodel.rxns),wrongNames));
for n = 1:length(rightRxns)
    if ismember(rightRxns(n),Rmodel.rxns)
        rxnNameID=strmatch(rightRxns(n),Rmodel.rxns,'exact');
        rightRxnEquations(n,1)=Rmodel.rxnEquations(rxnNameID);
        rightRxnNames(n,1)=Rmodel.rxnNames(rxnNameID);
        rightRxns(n,1)=Rmodel.rxns(rxnNameID);
        rightRev(n,1)=Rmodel.rev(rxnNameID);
        rightLB(n,1)=Rmodel.lb(rxnNameID);
        rightUB(n,1)=Rmodel.ub(rxnNameID);
        rightSubSystem(n,1)=Rmodel.subSystems(rxnNameID);
        rightGrRules(n,1)=Rmodel.grRules(rxnNameID);
    elseif ismember(rightRxns(n),Dmodel.rxns)
        rxnNameID=strmatch(rightRxns(n),Dmodel.rxns,'exact');
        rightRxnEquations(n,1)=Dmodel.rxnEquations(rxnNameID);
        rightRxnNames(n,1)=Dmodel.rxnNames(rxnNameID);
        rightRxns(n,1)=Dmodel.rxns(rxnNameID);
        rightRev(n,1)=Dmodel.rev(rxnNameID);
        rightLB(n,1)=Dmodel.lb(rxnNameID);
        rightUB(n,1)=Dmodel.ub(rxnNameID);
        rightSubSystem(n,1)=Dmodel.subSystems(rxnNameID);
        rightGrRules(n,1)=Dmodel.grRules(rxnNameID);
    end
end

%%% DUPLICATE REACTION HANDLING PROCEDURE
orthologsRemove = intersect(dictionary(:,2),Dmodel.genes);
tmpGPRs = Dmodel.grRules;
for n = 1:length(orthologsRemove)
    tmpGPRs = strrep(tmpGPRs,orthologsRemove(n),'False');
end
for n = 1:length(doubleRxns(:,1))
    rxnID1=strmatch(doubleRxns(n,1),rightRxns,'exact');
    rxnID2=strmatch(doubleRxns(n,2),Dmodel.rxns,'exact');
    tmpGPR = tmpGPRs{rxnID2};
    if ~isempty(tmpGPR)
        % brackets
        tmpGPR = regexprep(tmpGPR,'(.*)','\($1\)');
        tmpGPR = sprintf('''%s''',tmpGPR);
        [~,tmpGPR] = system(sprintf('python booleanSderenator.py %s',tmpGPR));
        % /brackets
        tmpGPR = regexprep(tmpGPR,'\((.*)\)','$1');
        rightGrRules(rxnID1)=strcat({'('},rightGrRules(rxnID1),{' or '},tmpGPR,{')'});
    end
end

[fused_model] = createModel(rightRxns,rightRxnNames,rightRxnEquations,rightRev,rightLB,rightUB,rightSubSystem,rightGrRules,'','');

%%RIPULIRE MODELLO: 
%unwantedGenes = setdiff(intersect(fused_model.genes,cat(1,recRemoveGenes,setdiff(donorModel.genes,addRepliconGenes))),'UnknownX');
%if ortFlag
%    badGenes = intersect(fused_model.genes,badGenes);
%    unwantedGenes = unique(cat(1,unwantedGenes,badGenes));
%end

%%Adjust the GPRs and stuff
%[fused_model,~,~,~] = deleteModelGenes(fused_model,unwantedGenes);

%[fused_model, wrongGPRs] = adjustFusedGPRs(fused_model);

%Remove safety 'X's
fused_model.grRules = strrep(fused_model.grRules,'X','');
fused_model.genes = strrep(fused_model.genes,'X','');
%wrongGPRs = strrep(wrongGPRs,'X','');
%wrongGPRs = reshape(wrongGPRs',[length(wrongGPRs)/2,2]);

%Check duplicate reactions
[~,removed]=checkDuplicateRxn(fused_model,2);
duplicateRxns=removed;
for n=1:length(removed)
    rxnID=strmatch(removed(n),fused_model.rxns);
    for a=1:length(fused_model.S)
        if rxnID~=a
            if isequal(fused_model.S(:,rxnID),fused_model.S(:,a))
                duplicateRxns(n,2) = fused_model.rxns(a);
                duplicateRxns(n,3) = strcat(fused_model.grRules(a),{' or '},fused_model.grRules(rxnID));
            end
        end
    end
end

swappedModel=fused_model;
end