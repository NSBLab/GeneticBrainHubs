filename = 'HCPparc20dens_4modelling.mat';

load(filename)

Afull = double(GrFA>0);

NNodesfull = length(Afull);
   
IndividualDists = zeros(NNodesfull,NNodesfull,length(individualCOORDs));

for j = 1:length(individualCOORDs)
   IndividualDists(:,:,j) = pdist2(individualCOORDs{j},individualCOORDs{j});
end
   
Dfull = mean(IndividualDists,3);

% We only run the model on one hemisphere


LHemi = NNodesfull/2;
A = Afull(1:LHemi,1:LHemi);
D = Dfull(1:LHemi,1:LHemi);

Aseed = zeros(length(A));

%% set parameters
% set limits for eta (geometric) and gamma (topological) parameters
eta = [-15,0];
gam = [-4,4];
% parameters related to the optimization
pow = 2;        % severity
nlvls = 5;      % number of steps
nreps = 2000;    % number of repetitions/samples per step
% specify model type and whether power-law or exponential (we really only
% used power-law).
modeltype = 'neighbors';
modelvar = [{'powerlaw'},{'powerlaw'}];
% number of edges
m = nnz(A)/2;

%% sample networks optimising for energy
% Function to use to calculate model fit
optimfunc = 'energy';

models = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

for i = 1:13
    
    modeltype = models{i};

    [E,K,N,P,C] = run_gen_model(A,Aseed,D,m,modeltype,modelvar,nreps,nlvl,etalim,gamlim,pow,optimfunc);
    
    [~,I] = min(E);
    
    NetsFromBestParams = generative_model(Aseed,Dist,m,modeltype,modelvar,repmat(P(I,:),100,1));
    
    [NetsFromBestParams_E,NetsFromBestParams_K,NetsFromBestParams_C] = fcn_eval_fits(A,Dist,NetsFromBestParams);   
    
    save(['Group_',modeltype,'_optenergy.mat'],'E','K','N','P','C',...
        'NetsFromBestParams','NetsFromBestParams_E','NetsFromBestParams_K','NetsFromBestParams_C','-v7.3')
    
end

%% sample networks optimising for correlation
% Function to use to calculate model fit
optimfunc = 'correlation';

models = {'sptl','neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod','com'};

for i = 1:13
    
    modeltype = models{i};

    [E,K,N,P,C] = run_gen_model(A,Aseed,D,m,modeltype,modelvar,nreps,nlvl,etalim,gamlim,pow,optimfunc);
    
    [~,I] = min(E);
    
    NetsFromBestParams = generative_model(Aseed,Dist,m,modeltype,modelvar,repmat(P(I,:),100,1));
    
    [NetsFromBestParams_E,NetsFromBestParams_K,NetsFromBestParams_C] = fcn_eval_fits(A,Dist,NetsFromBestParams);   
    
    save(['Group_',modeltype,'_optcorr.mat'],'E','K','N','P','C',...
        'NetsFromBestParams','NetsFromBestParams_E','NetsFromBestParams_K','NetsFromBestParams_C','-v7.3')
    
end