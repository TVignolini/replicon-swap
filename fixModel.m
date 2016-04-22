% generates a .rxnGeneMat and a .rules field for the input model

function fixedModel = fixModel(model)

tmpMat = sparse(length(model.grRules),length(model.genes));
tmpRules = regexprep(model.grRules,' and ',' \\ ');
tmpRules = regexprep(tmpRules,' or ',' \$ ');
for n = 1:length(model.grRules)
    tmpGenes = regexp(tmpRules(n),'[^\(\) \\\$]*','match');
    if ~isempty(tmpGenes{1})
        for m = 1:length(tmpGenes{1,:})
            tmpGene = tmpGenes{:};
            tmpGene = tmpGene{m};
            tmpID = strmatch(tmpGene,model.genes,'exact');
            tmpMat(n,tmpID) = 1;
            tmpID = num2str(tmpID);
            tmpRules(n) = regexprep(tmpRules(n),'(??@tmpGene)','x\(${tmpID}\)');
        end
    end
end

tmpRules = regexprep(tmpRules,' \\ ',' \& ');
tmpRules = regexprep(tmpRules,' \$ ',' \| ');

model = setfield(model,'rxnGeneMat',tmpMat);
model = setfield(model,'rules',tmpRules);

fixedModel = model;

end