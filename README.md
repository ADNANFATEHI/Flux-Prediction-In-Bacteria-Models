# Flux-Prediction-In-Bacteria-Models

This project compares flux distribution pattern on gene knockouts in Bacteria using six different Prediction Models.The Prediction models included are Flux Balance Analysis (FBA), Minimization of Metabolic Adjustments (MOMA),Metabolic Centric MOMA, Regulatory on/off Minimization (ROOM), Metabolic Centric ROOM and Relatch.

Requirements for using the project
1.	Install COBRA Toolbox written in Matlab. The link is here https://opencobra.github.io/
2.	An optimizer is to be installed to run these codes. The default solver with COBRA Toolbox is glpk (for LP and MILP).Solvers like    TOMLAB, IBM ILOG CPLEX, GUROBI, or MOSEK can be installed for faster simulation. The link is here https://github.com/opencobra/cobratoolbox/blob/master/docs/source/installation/solvers.md
3.	The RELATCH code was downloaded from http://reedlab.che.wisc.edu/
4.	Install heat map add on for Matlab from https://in.mathworks.com/matlabcentral/fileexchange/24253-customizable-heat-maps

Data Inputs required for New Analysis file are as follows

INPUT
1.	model (the following fields are required - others can be supplied)
    1) S-->Stoichiometric matrix
    2) b-->Right hand side = dx/dt
    3) c-->Objective coefficients
    4) lb-->Lower bounds
    5) ub-->Upper bounds

2.	MFA_Flux->Structure file containing data of Wild type Inner central Metabolism Reactions
    1) rxns-->Name of the Rxns
    2) val-->Flux values
    3) err-->Error in flux values
  
3.	External_Flux->Structure file containing data of Wild type External Rxns 
    1) rxns-->Name of the Rxns 
    2) val-->Flux values 
    3) err-->Error in flux values
  
4.	ReqRxnList->List of desired Reactions to be analysed

5.	geneList->List of genes to be knocked out

6.	expe->Experimental flux values of Mutant(Geneknockout) type

7.	expew->Experimental flux values of Wild type.

