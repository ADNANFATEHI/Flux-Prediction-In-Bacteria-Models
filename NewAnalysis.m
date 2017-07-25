%
% This script gives a comparison of different flux prediction models on geneknockout.
%INPUT
%1)model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%
%2)MFA_Flux->Structure file containing data of Wild type Inner central
%Metabolism Reactions
%   rxns        Name of the Rxns
%   val         Flux values
%   err         Error in flux values
%
%3)External_Flux->Structure file containing data of Wild type External
%Rxns
%   rxns        Name of the Rxns
%   val         Flux values
%   err         Error in flux values
%
%4)ReqRxnList->List of desired Reactions to be analysed
%
%5)geneList->List of genes to be knocked out
%
%6)expe->Experimental flux values of Mutant(Geneknockout) type
%
%7)expew->Experimental flux values of Wild type.
%
%OUTPUT
%1)Fluxes_FBA->Flux values predicted by FBA model
%2)MOMA->Flux values predicted by MOMA model
%3)MOMAME->Flux values predicted Metabolic Centric MOMA model
%4)Flux->Flux values predicted by ROOM model
%5)RELATCH->Flux values predicted by RELATCH model
%6)ROOMME->Flux values predicted by Metabolic Centric ROOM model
%7)Provides Heat Maps fro the above models
%8)Provides Average Error of the different models in form of a bar graph
%
%Created by Adnan Fatehi
%Date 21/07/17
%
global wildf;
load 'RELATCH_iAF1260_Example.mat' %Upload model here
nobj=findRxnIDs(model,{'Ec_biomass_iAF1260_core_59p81M'}); %Define objective function
% model.ub(nobj)=0.1;   %For Continuous Culture set up upper bounds
% model.lb(nobj)=0.1;   %For Continuous Culture set up lower bounds

load MFA_Flux;          %load MFA file here
load External_Flux;     %load External Flux here
load 'ReqRxnList';      %load required ReqRxnList here
load 'geneList.mat';    %load geneList here    
load 'expe.mat';        %load experiment mutant flux values here
load 'expew.mat';       %load experiment wild flux values here 

%Updating Optimize tools here 
changeCobraSolver('ibm_cplex','MILP',1);
changeCobraSolver('ibm_cplex','LP',1);
changeCobraSolver('ibm_cplex','QP',1);

%Solver to find Reference solution
solver='cplex_direct';
modelRef = model;
solutionRef = RELATCH_Reference(modelRef,Gene_Expression,External_Flux,MFA_Flux,solver);
wildf=solutionRef.w;

%Updating Optimize tools here 
changeCobraSolver('ibm_cplex','MILP',1);
changeCobraSolver('ibm_cplex','LP',1);
changeCobraSolver('ibm_cplex','QP',1);
    
ng=length(geneList);        %counts the number of genes to be deleted
nRxns=length(model.rxns);   %counts the number of reactions


% model.ub(nobj)=0.22;
% model.lb(nobj)=0.22;


for i=1:ng
    
    %FBA STARTS
    
    tempgene=geneList(i);
    [all_Fluxes_FBA(:,i) fval(i) ] = FBA1( model,tempgene,nobj);
    
    disp('FBA DONE');
    %FBA ENDS
    
    %MOMA STARTS
        tempgene=geneList(i);
        [all_Fluxes_MOMA(:,i) fval] = MOMA1( model,tempgene);
        
         disp('MOMA DONE');
    %MOMA ENDS
   
        
    %RELATCH STARTS
    alpha=10;
    gamma=1.1;
    downRegFraction=0;
    tempgenerel=geneList(i);
    [tempmodelrel,hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(model,tempgenerel,downRegFraction);

    solver='cplex_direct';
    try
       solutionPert {i}= RELATCH_Perturbed(tempmodelrel,solutionRef,alpha,gamma,solver);
       all_Fluxes_RELATCH(:,i)=solutionPert{i}.v;         
    catch
        solnPert(i)=1;
    end
    disp('RELATCH DONE');
    %RELATCH ENDS

    %MOMAME STARTS
    
    tempgene1=geneList(i);
    [all_Fluxes_MOMAME(:,i)] = MOMAME(model,tempgene1);
    
     disp('MOMAME DONE');
    %MOMAME ENDS

    %ROOM STARTS
    tempgene=geneList(i);
     try
    
     [all_Fluxes_ROOM(:,i) fval(i)] = ROOM( model,tempgene);
   
     catch
     all_Fluxes_ROOM(:,i)=0;
     end
    disp('ROOM DONE');
    %ROOMME STARTS
    try
    tempgene1=geneList(i);
    [all_Fluxes_ROOMME(:,i)] = ROOMME1(model,tempgene1);
    catch
        all_Fluxes_ROOMME(:,i)=0;
    end
   
    %ROOMME ENDS
    
disp(i);
end
 
goodRxns=intersect(model.rxns,ReqRxnList);  %selecting goodRxns from model
ngood=length(goodRxns);                     %find the numbe of good Reactions
idgood=findRxnIDs(model,goodRxns);          %finding ids of good Reactions
wilds=wildf(idgood);

%Extracting the required reactions the model and generate heat map
Fluxes_FBA=all_Fluxes_FBA(idgood,:);       
clabel = arrayfun(@(x){sprintf('%0.1f',x)},Fluxes_FBA);
heatmap(Fluxes_FBA,geneList,ReqRxnList,clabel,'TickAngle',90,'ShowAllTicks', true,'FontSize',10, 'TextColor', 'black',...
  'UseLogColormap', true,'GridLines', ':');
title('FBA');
figure;

MOMA=all_Fluxes_MOMA(idgood,:);
clabel = arrayfun(@(x){sprintf('%0.1f',x)},MOMA);
heatmap(MOMA,geneList,ReqRxnList,clabel,'TickAngle',90,'ShowAllTicks', true,'FontSize',10, 'TextColor', 'black',...
  'UseLogColormap', true,'GridLines', ':');
title('MOMA');
figure;

MOMAME=all_Fluxes_MOMAME(idgood,:);
clabel = arrayfun(@(x){sprintf('%0.1f',x)},MOMAME);
heatmap(MOMAME,geneList,ReqRxnList,clabel,'TickAngle',90,'ShowAllTicks', true,'FontSize',10, 'TextColor', 'black',...
  'UseLogColormap', true,'GridLines', ':');
title('MCMOMA');
figure;

Flux= all_Fluxes_ROOM(idgood,:);
clabel = arrayfun(@(x){sprintf('%0.1f',x)},Flux);
heatmap(Flux,geneList,ReqRxnList,clabel,'TickAngle',90,'ShowAllTicks', true,'FontSize',10, 'TextColor', 'black',...
  'UseLogColormap', true,'GridLines', ':');
title('ROOM');
figure;

RELATCH=all_Fluxes_RELATCH(idgood,:);
clabel = arrayfun(@(x){sprintf('%0.1f',x)},RELATCH);
heatmap(RELATCH,geneList,ReqRxnList,clabel,'TickAngle',90,'ShowAllTicks', true,'FontSize',10, 'TextColor', 'black',...
  'UseLogColormap', true,'GridLines', ':');
title('RELATCH');
figure;

ROOMME=all_Fluxes_ROOMME(idgood,:);
clabel = arrayfun(@(x){sprintf('%0.1f',x)},ROOMME);
heatmap(ROOMME,geneList,ReqRxnList,clabel,'TickAngle',90,'ShowAllTicks', true,'FontSize',10, 'TextColor', 'black',...
  'UseLogColormap', true,'GridLines', ':');
title('ROOMME');
figure;

%Convert cells to double
expe=cell2mat(expe);
expew=cell2mat(expew);

%Finding Absolute Eroor 
FBAerror=((Fluxes_FBA-expe).^2);
nFBA=(FBAerror~=0);
sumFBA=sum(FBAerror)./sum(nFBA);

MOMAerror=((MOMA-expe).^2);
nMOMA=(MOMAerror~=0);
sumMOMA=sum(MOMAerror)./sum(nMOMA);

MOMAMEerror=((MOMAME-expe).^2);
nMOMAME=(MOMAMEerror~=0);
sumMOMAME=sum(MOMAMEerror)./sum(nMOMAME);


ROOMerror=((Flux-expe).^2);
nROOM=(ROOMerror~=0);
sumROOM=sum(ROOMerror)./sum(nROOM);

ROOMMEerror=((ROOMME-expe).^2);
nROOMME=(ROOMMEerror~=0);
sumROOMME=sum(ROOMMEerror)./sum(nROOMME);

RELATCHerror=((RELATCH-expe).^2);
nRELATCH=(RELATCHerror~=0);
sumRELATCH=sum(RELATCHerror)./sum(nRELATCH);

wilderror=((wilds-expew).^2);
nwild=(wilderror~=0);
sumwild=sum(wilderror)./sum(nwild);

%Displaying Bar Graphs of errors
ax = gca;
ax.XTickLabels = {'FBA','MOMA','MOMA MET','ROOM','ROOM MET','RELATCH';};
ax.XTickLabelRotation = 45;
y = [sumFBA sumMOMA sumMOMAME sumROOM sumROOMME sumRELATCH];    
b=bar(ax,y); 
