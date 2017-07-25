function [ all_Fluxes_FBA fvalFBA ] = FBA1( model ,tempgene, nobj)
%
%FBA1 makes use of linear programing to perform Flux Balance Analysis.
%
%
%   [ x fval ] = FBA1( model ,tempgene)
%
%   Mathematic Model
%   Maximize Z=C*V
%   Subject to Constraints	S.V=0
%		lb<V<ub
%   where lb? Lower bound of flux
%         ub? Upper bound of flux
%	       V? Flux value in mmol gDW–1 h–1
%
%   For biomass production, the value of c=1 and for remaining reactions c=0 in objective function. 
%
%INPUTS
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds 
%
% tempgene       Gene to be deleted.
% nobj           RxnID of the objective function
%OUTPUTS
% fvalFBA           Objective function value.
% x              Value of Flux passing through each reaction.
%
%Notes
%1)Used removeInfeasibleLoops function file to removes high flux values.
%  Function file uploaded.
%
%Created by Adnan Fatehi
%Date 21/07/17

nRxns=length(model.rxns); %Finds the number of Rxns in model.
%Checks weather gene is present or not. If no gene present perfoms FBA for wild type strain
n=isempty(tempgene);
options=cplexoptimset;
options.Algorithms='Interior-point';
options.Display='off';

if(n==1)
    A=(model.S);                                        %Represents the S matrix
    B=model.b;                                          %Represents the RHS of the eqn.
    f=zeros(nRxns,1); 
    f(nobj)=1;
    [y, fvalFBA]=cplexlp(-f,[],[],A,B,model.lb,model.ub,[],options);  %Liner programing S.v=0
    [all_Fluxes_FBA,Corrected]=removeInfeasibleLoops(model,y);       %Remove infeasible Loops
    all_Fluxes_FBA(abs(all_Fluxes_FBA)>100)=0;
% Flux Calculation when gene deletion is there    
else
    tempmodel=model;                                        %Temporary model created
    downRegFraction=0;

% Model is stored in tempmodel after deleting gene. Gene is deleted using deleteModelGene fuction file.

 [tempmodel,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(tempmodel,tempgene,downRegFraction);
 A=(tempmodel.S);                                       %Represents the S matrix
 B=tempmodel.b;                                             %Represents the RHS of the eqn.
 f=zeros(nRxns,1);
 f(nobj)=1;
 [y, fvalFBA]=cplexlp(-f,[],[],A,B,tempmodel.lb,tempmodel.ub,[],options);     %Liner programing S.v=0
 [all_Fluxes_FBA,Corrected]=removeInfeasibleLoops(model,y);  %Remove infeasible Loops
 all_Fluxes_FBA(abs(all_Fluxes_FBA)>100)=0;
end


end


