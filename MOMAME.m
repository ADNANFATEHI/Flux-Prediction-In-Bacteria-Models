function [x] = MOMAME(model,tempgene1)
%
%MOMAME uses Metabolic Centric Approach to perform Minimization of
%Metabollic Adjustments (MOMA). It makes use of Linear Least squares to
%minimize the Euclidean distance between metabollites of wild and mutant
%flux.
%
%
%[x] = MOMAME( model ,tempgene1)
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
%   Considering Summation of Incoming flux for each metabolite
%           ?=0.5 ?s?v
%
%   Minimize Z= ?(CX-D)^2  
%   C represents stochiometric matrix containing (+) ve coefficients of metabolites
%   D represents matrix of wild type metabolite flux.
%
%   Subject to Constraints	S.V=0
%		lb<V<ub
%   where lb? Lower bound of flux
%         ub? Upper bound of flux
%         wi?Wild Type Flux in mmol gDW–1 h–1
%          V? Flux value in mmol gDW–1 h–1
%
%INPUTS
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds 
%
% tempgene1       Gene to be deleted.
%
%
%OUTPUTS
% x              Value of Flux passing through each reaction.
%
%Notes
%1)In this wild type flux is calculated using FBA1 and stored in wildf. It is declared as a global variable.
%2)Used convertToIrreversible function from CobraToolbox to convert model to
%  irreversible format.
%3)Used convertRevFluxDistribution function file to convert reversible wild
%  type flux calcuclated using FBA to irreversible flux. Function file
%  uploaded.
%4)Used convertIrrevFluxDistribution function file from CobraToolbox to
%  convert irreversible flux to reversible flux.
%5)Used removeInfeasibleLoops function file to removes high flux values.
%  Function file uploaded.
%
%Created by Adnan Fatehi
%Date 21/07/17

global wildf
modelRef=model;
nr=length(model.rxns);

%Converting reversible model to irreversible model
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(modelRef);

nrxns=length(modelIrrev.rxns);
nmet=length(modelIrrev.mets);

Flux_FBAIrrev = convertRevFluxDistribution(wildf,rev2irrev,nr);%Converting reversible flux to irreversible

CW=modelIrrev.S;

%Considering Positive Fluxes
for i=1:nmet
    for j=1:nrxns
      if (modelIrrev.S(i,j)<0)
         CW(i,j)=0;
      end
    end
end

wild=(CW)*(Flux_FBAIrrev); 

%Knocking out gene
tempmodel=modelIrrev;
downRegFraction=0;
[tempmodel,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(tempmodel,tempgene1,downRegFraction);

%Scaling down the matrices

smallestMet=1e-3;
zeroRows=find(wild<smallestMet);
nonzeroRows=setdiff(1:nmet,zeroRows);
CW(zeroRows,:)=0;
CW(nonzeroRows,:)=CW(nonzeroRows,:)./repmat(wild(nonzeroRows),1,nrxns);
d=ones(nmet,1);
d(zeroRows)=0;

%Setting up equality Constraint AX=B
Aeq=(tempmodel.S);             
Beq=tempmodel.b;

%Represents the positive coefficients of S matrix
C=CW;

options=optimoptions('lsqlin','Algorithm','interior-point','MaxIter',200);
x= lsqlin(C,d,[],[],Aeq,Beq,tempmodel.lb,tempmodel.ub,[],options);%CX-D square

x= convertIrrevFluxDistribution(x,matchRev); %Converts Irreversible FLux distribution to reversible
[x,Corrected]=removeInfeasibleLoops(model,x); %Remove Infeasible Loops

end

