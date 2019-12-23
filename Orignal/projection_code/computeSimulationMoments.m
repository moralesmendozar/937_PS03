function [simseries, varnames] = computeSimulationMoments(simseries, varnames)

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:length(varnames)
    if isempty(indexmap.get(varnames{i}))
        indexmap.put(varnames{i},i);
    end
end

%--------------------------------------------------------------------------
% transformations of raw model output
%--------------------------------------------------------------------------

% list of indices
loglist=Helpers.makeListFromNames(indexmap,{'R','q','cF','cG','kGpol'});
multlist=Helpers.makeListFromNames(indexmap,{'muF','lamF','lamG'});

% conversion of log-values
simseries(:,loglist)=exp(simseries(:,loglist));
% conversion of multipliers
simseries(:,multlist)=max(simseries(:,multlist),0).^(1/3);

% ---------------------------------------------------------------------
% borrower
% ---------------------------------------------------------------------
bind_muF= simseries(:,indexmap.get('muF'))>0;
q= simseries(:,indexmap.get('q'));
wF= simseries(:,indexmap.get('wF'));
wFsh= wF./q;
eff_theta= simseries(:,indexmap.get('min_q'))./q;

simseries=[simseries(:,2:end),bind_muF,wFsh,eff_theta];

varnames_add={'bind_muF','wFsh','eff_theta'};
varnames=[varnames(2:end), varnames_add];

end



