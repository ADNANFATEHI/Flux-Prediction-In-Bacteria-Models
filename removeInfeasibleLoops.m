function [newFluxes,Corrected]=removeInfeasibleLoops(model,oldFluxes)
% remove thermodynamically infeasible loops from flux distributions
%set upper and lower bounds
for i=1:length(model.rxns)
    if (strcmp(model.rxns(i),'Ec_biomass_iAF1260_core_59p81M'))
%         strcmp(model.subSystems(i),'Exchange') ||
        model.lb(i)=oldFluxes(i);
        model.ub(i)=oldFluxes(i);
    elseif oldFluxes(i)>0
        model.lb(i)=0;
        model.ub(i)=oldFluxes(i);
    else
        model.ub(i)=0;
        model.lb(i)=oldFluxes(i);
    end
    
end
% minimize the absolute sum of the fluxes
sol = optimizeCbModel(model,[],'one');
if ~isempty(sol.x)
    newFluxes=sol.x;
    Corrected=1;
else
    newFluxes=oldFluxes;
    Corrected=0;
end
end

