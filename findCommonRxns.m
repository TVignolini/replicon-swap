%%% Finds common reactions between 2 models %%%
%%% ancora problemi per quanto riguarda reazioni che hanno stessi
%%% metaboliti ma coefficienti diversi. Le calcola come uguali. Questo per
%%% il problema dell'approssimazione di cui prima

function [commonRxns] = findCommonRxns(model1,model2)

%Create a single big S matrix in which to compare the reactions
Smets = unique(cat(1,model1.mets,model2.mets));
Srxns = cat(1,model1.rxns,model2.rxns);
Stmp = sparse(length(Smets),length(Srxns));
for n = 1:length(model1.rxns)
    rxnID = strmatch(Srxns(n),model1.rxns,'exact');
    Scolumn = model1.S(:,rxnID);
    for m = 1:length(Scolumn)
        if any(Scolumn(m))
            metID = strmatch(model1.mets(m),Smets,'exact');
            Stmp(metID,n) = Scolumn(m);
        end
    end
end
for n = (length(model1.rxns)+1):length(Srxns)
    rxnID = strmatch(Srxns(n),model2.rxns,'exact');
    Scolumn = model2.S(:,rxnID);
    for m = 1:length(Scolumn)
        if any(Scolumn(m))
            metID = strmatch(model2.mets(m),Smets,'exact');
            Stmp(metID,n) = Scolumn(m);
        end
    end
end

%Compare matrix columns to determine equivalent reactions
commonRxns = cell(length(Srxns),2);
for n = 1:(length(Srxns)-1)
    for m = n+1:length(Srxns)
        if isequal(Stmp(:,n),Stmp(:,m))
            commonRxns(n,1) = Srxns(n);
            commonRxns(n,2) = Srxns(m);            
        end
    end
end
commonRxns = commonRxns(~cellfun('isempty',commonRxns));
commonRxns = reshape(commonRxns,length(commonRxns)/2,2);
end