function [ x fval] = MOMA1( model,tempgene )
%
%MOMA1 makes use of quadratic programing to perform Minimization of
%Metabollic Adjustments (MOMA) to minimize the Euclidean Distance between
%Wild type and Mutant Flux.
%
%
%[ x fval ] = MOMA1( model ,tempgene)
%
%   Mathematic Model,
%   Minimize Z= ?(wi-vi)2
%   Z=vi^2-2*vi*wi	(wi2 neglected as its constant and does not affect the minimization)
%   Q represents coefficient of vi2
%   L represents coefficient of vi that is 2wi
%   
%   Subject to Constraints	S.V=0
%                           lb<V<ub
%   where lb? Lower bound of flux
%         ub? Upper bound of flux
%         wi?Wild Type Flux in mmol gDW–1 h–1
%         V? Flux value in mmol gDW–1 h–1
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
% fval           Distance between wild type and mutant type flux.
% x              Value of Flux passing through each reaction.
%
%Notes
%1)In this wild type flux is calculated using FBA1 and stored in wildf. It is declared as a global variable.
%2)Used removeInfeasibleLoops function file to removes high flux values.
%  Function file uploaded.
%
%Created by Adnan Fatehi
%Date 21/07/17
global wildf

 nRxns=length(model.rxns); %Finds the number of Rxns in model.
 tempmodel=model;
 downRegFraction=0;
 %Deleting a gene
 [tempmodel,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(tempmodel,tempgene,downRegFraction);
 %Setting up equality Constraint AX=B
 A=(tempmodel.S);
 B=tempmodel.b;
 
 const=-2*wildf;
 f=zeros(nRxns,1); %Setting Objective Function
 f(:,1)=const;
 H=2*eye(nRxns);
 [x, fval]=quadprog(H,f,[],[],A,B,tempmodel.lb,tempmodel.ub); %Solving quadratic prog (w-v)^2 
 [x,Corrected]=removeInfeasibleLoops(tempmodel,x);  %remove infeasible loops
  
end

