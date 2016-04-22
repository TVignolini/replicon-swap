% Applies BooleanSderenator.py

function sderenatedModel = sderenateModel(model)

% The more brackets the better
model.grRules = regexprep(model.grRules,'(.*)','\($1\)');

for n = 1:length(model.grRules)
    if ~isempty(model.grRules{n})
        newRule = sprintf('''%s''',model.grRules{n});
        [~,newRule] = system(sprintf('python booleanSderenator.py %s',newRule));
        model.grRules(n) = cellstr(newRule);
    end
end

% Too many brackets
model.grRules = regexprep(model.grRules,'\((.*)\)','$1');

sderenatedModel = model;

end