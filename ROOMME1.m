function [ all_Fluxes_ROOMME3] =ROOMME1( model, tempgene1 )
%
%ROOMME1 uses Metabolic Centric Approach to perform Regulatory on/off minimization.
%It makes use of mix integer linear programing to minimize the Hamming
%distance between metabolites of wild and mutant type flux.
%
%
%[ all_Fluxes_ROOMME3 ] = ROOMME1( model ,tempgene1)
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
%   Mathematic Model
%   Minimize ?yi (i=1 to N)
%   Subject to Constraints
%   S?v=0   v(min)<v<v(max)
%       for 1<i<n
%   ?vi- yi*(?v min,i - wild(i,l)) >wild(l)
%   ?vi- yi*(?v max,i - wild(i,u)) < wild(u)
%        yi?(0,1)
%   wiu  = ?wi +? |?wi|+?
%   wil  = ?wi -? |?wi|-?
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
%
%
%OUTPUTS
% all_Fluxes_ROOMME3    Value of Flux passing through each reaction.
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

global wildf;
modelRef=model;
nr=length(model.rxns);
%Converting reversible model to irreversible model
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(modelRef);

nrxns=length(modelIrrev.rxns);
nmet=length(modelIrrev.mets);


Flux_FBAIrrev = convertRevFluxDistribution(wildf,rev2irrev,nr); %Converting reversible flux to irreversible

CW=modelIrrev.S;
%Considering Positive Fluxes
for i=1:nmet
    for j=1:nrxns
      if (modelIrrev.S(i,j)<0)
         CW(i,j)=0;
      end
    end
end
wildm=(CW)*(Flux_FBAIrrev);

    
eta=0.01;
delta=0.1;
wildmu=wildm+delta*abs(wildm)+eta;     %constant RHS uperbound terms 
wildml=wildm-delta*abs(wildm)-eta;     %constant RHS lower bound term
cmax=CW*modelIrrev.ub-wildmu;          %constant LHS uperbound terms 
cmin=CW*modelIrrev.lb-wildml;          %constant LHS lower bound term


%Knocking out gene
tempmodel=modelIrrev;
downRegFraction=0;
[tempmodel,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(tempmodel,tempgene1,downRegFraction);

%Setting up equality Constraint AX=B
Aeq=(tempmodel.S);             %Represents the S matrix
Beq=tempmodel.b;

f=[zeros(nrxns,1);ones(nmet,1)];  %Objective Function  
intcon=(nrxns+1):(nrxns+nmet);    %Integer values in the objective function

A8=(zeros(nmet,nmet));           % create new integer constraints

% add new integer constraints
Aeq2=[Aeq A8];
   
%defing constraint1
Aineq1=[CW  -diag(cmax)];       
bineq1=wildmu;

%defing constraint2    
Aineq2=[-CW diag(cmin)];  
bineq2=-wildml;
    
%combining einequality constraints   
Aineq=[Aineq1;Aineq2];    
bineq=[bineq1;bineq2];
       
%add new integer constraints to the bounds 
    
y=[tempmodel.lb;zeros(nmet,1)];    
z=[tempmodel.ub;ones(nmet,1)];

[all_Fluxes_ROOMME]= intlinprog(f,intcon,Aineq,bineq,Aeq2,Beq,y,z); %MILP program solver
all_Fluxes_ROOMME1=all_Fluxes_ROOMME(1:nrxns);  %selcts the Rxns and not the constraints     
all_Fluxes_ROOMME2= convertIrrevFluxDistribution(all_Fluxes_ROOMME1,matchRev);  %Converts teh selected rxns into reversible
[all_Fluxes_ROOMME3,Corrected]=removeInfeasibleLoops(model,all_Fluxes_ROOMME2); %Remove Infeasible Loops
    
end


