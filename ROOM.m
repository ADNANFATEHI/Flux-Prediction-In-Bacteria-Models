function [all_Fluxes_ROOM fval] = ROOM(model,tempgene)
%
%ROOM makes use of mix integer linear programing to perform Regulatory on/
%off minimiztion.It minimizes the number of significant flux changes in the
%model by minimizing the hamming distance between wild and mutant type
%flux.
%
%
%[all_Fluxes_ROOM] = ROOM( model ,tempgene)
%
%   Mathematic Model
%   Minimize ?yi (i=1 to N)
%   Subject to Constraints
%   S?v=0   v(min)<v<v(max)
%       for 1<i<n
%   vi- yi*(v min,i - wild(i,l)) >wild(l)
%   vi- yi*(v max,i - wild(i,u)) < wild(u)
%        yi?(0,1)
%   wiu  = wi +? |wi|+?
%   wil  = wi -? |wi|-?
%
%INPUTS
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds 
%
%tempgene       Gene to be deleted.
%
%OUTPUTS
%fval               Objective function value.
%all_Fluxes_ROOM    Value of Flux passing through each reaction.
%
%Notes
%1)In this wild type flux is calculated using FBA1 and stored in wildf. It is declared as a global variable.
%2)Used removeInfeasibleLoops function file to removes high flux values.
%  Function file uploaded.
%
%Created by Adnan Fatehi
%Date 07/21/17.

global wildf;
nRxns=length(model.rxns);       %Finds the number of Rxns in Model
nMets=length(model.mets);       %Finds the number of metabolites in model
tempmodel=model;                
downRegFraction=0;
% Model is stored in tempmodel after deleting gene. Gene is deleted using deleteModelGene fuction file.
[tempmodel,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(tempmodel,tempgene,downRegFraction);

%Setting up equality Constraint AX=B
Aeq=(tempmodel.S);             %Represents the S matrix
Beq=tempmodel.b;               %Represents the RHS of the eqn.

%calculation of constant
eta=0.01;
delta=0.1;
wildfu=wildf+delta*abs(wildf)+eta;      %constant RHS uperbound terms 
wildfl=wildf-delta*abs(wildf)-eta;      %constant RHS lowerbound terms 
cmax=tempmodel.ub-wildfu;               %constant LHS uperbound terms
cmin=tempmodel.lb-wildfl;               %constant LHS lowerbound terms



%Setting up objective function
f=[zeros(nRxns,1);ones(nRxns,1)];        %Objective Function
intcon=(nRxns+1):nRxns*2;                %Integer values in the objective function
A8=(zeros(nMets,nRxns));                 % create new integer constraints


% add new integer constraints
Aeq2=[Aeq A8];

%defing constraint1
Aineq1=[eye(nRxns) -diag(cmax)];
bineq1=wildfu;

%defing constraint1
Aineq2=[-eye(nRxns) diag(cmin)];
bineq2=-wildfl;

%combining inequality constraints
Aineq=[Aineq1;Aineq2];
bineq=[bineq1;bineq2];
%add new integer constraints to bounds

y=[tempmodel.lb;zeros(nRxns,1)];
z=[tempmodel.ub;ones(nRxns,1)];

[all_Fluxes_ROOM1 fval]= intlinprog(f,intcon,Aineq,bineq,Aeq2,Beq,y,z); %MILP program solver
all_Fluxes_ROOM2=all_Fluxes_ROOM1(1:nRxns);

[all_Fluxes_ROOM,Corrected]=removeInfeasibleLoops(model,all_Fluxes_ROOM2); %Remove Infeasible Loops

end
