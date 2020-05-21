function mask = makeRFPmask(nodeDeg, khub, Adj, j)

% j=1 gives mask or RICH links
% j=2 gives a mask for feeder links
% j=3 guves a mask for peripheral links
isHub = nodeDeg>khub; 
mask = nan(length(nodeDeg), length(nodeDeg)); 
        switch j
            case 1 % 'rich'
                mask(isHub,isHub) = 1;
            case 2 % 'feeder'
                mask(~isHub,isHub) = 1;
                mask(isHub,~isHub) = 1;
            case 3 % 'local'
                mask(~isHub,~isHub) = 1;
        end

        mask = mask.*Adj; 
        mask(mask==0) = NaN; 
end