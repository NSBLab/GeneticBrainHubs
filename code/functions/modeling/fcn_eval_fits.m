function [E,KS,C] = fcn_eval_fits(A,Dist,b)
n = length(Dist);
ndraw = size(b,2);
ka = sum(A,2);
ba = betweenness_bin(A)';
da = Dist(triu(A,1) > 0);
ca = clustering_coef_bu(A);
KS = zeros(ndraw,4);
C = zeros(ndraw,1);
for idraw = 1:ndraw
    a = zeros(n);
    a(b(:,idraw)) = 1;
    a = a + a';
    x = sum(a,2);
    y = betweenness_bin(a)';
    z = Dist(triu(a,1) > 0);
    zz = clustering_coef_bu(a);
    KS(idraw,1) = fcn_ks(x,ka);
    KS(idraw,2) = fcn_ks(y,ba);
    KS(idraw,3) = fcn_ks(z,da);
    KS(idraw,4) = fcn_ks(ca,zz);
    
    C(idraw) = corr(x,ka,'Type','Spearman');
end
E = max(KS,[],2);