function [W_thr,Wmeasure] = threshold_consistency_AA(Ws, p, CVmeasure, consThr)
%THRESHOLD_CONSISTENCY    Threshold edges ranked by consistency
%
%   W_thr = threshold_consistency(Ws, p);
%
%   This function "thresholds" the group mean connectivity matrix of a set
%   of connectivity matrices by preserving a proportion p (0<p<1) of the
%   edges with the smallest coefficient of variation across the group. All
%   other weights, and all weights on the main diagonal (self-self
%   connections) are set to 0.
%
%   Inputs: Ws,     N-by-N-by-M group of M weighted connectivity matrices
%           p,      proportion of weights to preserve
%                       range:  p=1 (all weights preserved) to
%                               p=0 (no weights removed)
%
%   Output: W_thr,  thresholded group mean connectivity matrix
%
% first, select only edges that are present in at least 50% subjects
LWs = logical(Ws);
consistWs = mean(LWs, 3);
edgesKeep = consistWs>=consThr;

Wmeasure = zeros(size(Ws,1), size(Ws,1));
if strcmp(CVmeasure, 'strength')
    Wmeasure(~edgesKeep) = -inf;
else
    Wmeasure(~edgesKeep) = inf;
end

for i=1:length(Wmeasure)
    for j=i+1:length(Wmeasure)
        
        % if the edge exists in 50% subjects consider only nonzero values
        % to calculate the measure of consistency
        if edgesKeep(i,j)
            edge = squeeze(Ws(i,j,:));
            %edge = squeeze(Ws(i,j,:));
            switch CVmeasure
                case 'RMD' % relative mean absolute difference
                    Wmeasure(i,j) = mad(edge,0)./mean(edge);
                case 'QCD' % quantile coefficient of dispersion
                    Wmeasure(i,j) = (quantile(edge, 0.75) - quantile(edge, 0.25))/(quantile(edge, 0.75) + quantile(edge, 0.25));
                case 'MAD' % median absolute deviation/median
                    Wmeasure(i,j) = mad(edge,1)./median(edge);
                    %                 case 'MADpoisson' % median absolute deviation/median
                    %
                    %                     L = var(edge)/100000; % divide by a big number, so the
                    %                     fL = floor(L);
                    %                     Wmeasure(i,j) = 2*exp(-L*(L^(fL+1))/(factorial(fL)));
                    %                 case 'MADpoisson0' % median absolute deviation/median
                    %                     edge = edge./1000;
                    %                     L = ((var(edge)^2 + mean(edge)^2)/mean(edge)) - 1;  % divide by a big number, so the
                    %                     fL = floor(L);
                    %                     Wmeasure(i,j) = 2*exp(-L*(L^(fL+1))/(factorial(fL)));
                case 'CV'
                    Wmeasure(i,j)=std(edge)./mean(edge); % coefficient of variation across the group
                case 'strength'
                    Wmeasure(i,j)= median(edge); % coefficient of variation across the group
                case 'CVmod' % quantile coefficient of dispersion
                    Wmeasure(i,j) = (quantile(edge, 0.75) - quantile(edge, 0.25))/quantile(edge, 0.5);
            end
        end
    end
end

% get the mean of existing edges
if strcmp(CVmeasure, 'strength')
    Ws(isnan(Ws)) = -inf;
else
    Ws(isnan(Ws)) = inf;
end
Wmean = nanmean(Ws,3);

Wmeasure = Wmeasure+Wmeasure';
if strcmp(CVmeasure, 'strength')
    Wmeasure(isnan(Wmeasure)) = -inf;
    Wmeasure(Wmeasure==0) = -inf;
else
    Wmeasure(isnan(Wmeasure)) = inf;
    Wmeasure(Wmeasure==0) = inf;
end

if strcmp(CVmeasure, 'strength')
    W_thr=threshold_arbmeasure(Wmean,Wmeasure,p); %
else
    W_thr=threshold_arbmeasure(Wmean,-Wmeasure,p); % th
end
end


