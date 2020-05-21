function [energyout,ksout,netsout,ptsout,corrout] = ...
    run_gen_model(A,Aseed,Dist,m,modeltype,modelvar,ndraw,nlvl,etalim,gamlim,pow,optimfunc)

indseed = find(triu(Aseed,1));
mseed = nnz(Aseed)/2;
mrem = m - mseed;
totalsamples = ndraw*nlvl;
ptsout = zeros(totalsamples,2);
energyout = zeros(totalsamples,1);
netsout = zeros(mrem,totalsamples);
ksout = zeros(totalsamples,4);

corrout = zeros(totalsamples,1);

eta = unifrnd(etalim(1),etalim(2),ndraw,1);
gam = unifrnd(gamlim(1),gamlim(2),ndraw,1);
ptsout(1:ndraw,:) = [eta,gam];

fprintf('level %i of %i\n',1,nlvl);

netstemp = generative_model(Aseed,Dist,m,modeltype,modelvar,ptsout(1:ndraw,:));

[energyout(1:ndraw),ksout(1:ndraw,:),corrout(1:ndraw)] = fcn_eval_fits(A,Dist,netstemp);

for i = 1:mseed
    netstemp(netstemp == indseed(i)) = [];
end
netstemp = reshape(netstemp,[mrem,ndraw]);
netsout(:,1:ndraw) = netstemp;

powvals = linspace(0,pow,nlvl);

for ilvl = 2:nlvl
    
    pow = powvals(ilvl);
    
    fprintf('level %i of %i\n',ilvl,nlvl);
    ind = 1:ndraw*(ilvl - 1);
    
    switch optimfunc
        case 'energy'
            fitfunc = energyout(ind);       
        case 'correlation'
            fitfunc = 1+(corrout(ind).*-1);
    end
    
    if strcmp(modeltype,'sptl')
        ptsnew = fcn_voronoi_select_sptl(ptsout(ind,1),fitfunc(ind),ndraw,etalim,pow);
    else
        ptsnew = fcn_voronoi_select(ptsout(ind,:),fitfunc(ind),ndraw,etalim,gamlim,pow);
    end
    
    indnew = (1 + (ilvl - 1)*ndraw):(ilvl*ndraw);
    
    netstemp = generative_model(Aseed,Dist,m,modeltype,modelvar,ptsout(1:ndraw,:));

    [energyout(indnew),ksout(indnew,:),corrout(indnew)] = fcn_eval_fits(A,Dist,netstemp);
    
    if strcmp(modeltype,'sptl')
        ptsout(indnew,1) = ptsnew;
    else
        ptsout(indnew,:) = ptsnew;
    end
    for i = 1:mseed
        netstemp(netstemp == indseed(i)) = [];
    end
    netstemp = reshape(netstemp,[mrem,ndraw]);
    netsout(:,indnew) = netstemp;
        
end