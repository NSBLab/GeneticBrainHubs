function [FitCurve,Residuals] = getFitCurve(Fit,distExpVect,c)

Dvect = distExpVect(:,1);
if nargin < 4
    xThresholds = 25;
end

switch Fit{1}
    
    case 'linear'
        FitCurve = c.p1*Dvect + c.p2;
    case 'exp'
        FitCurve = c.A*exp(-c.n*Dvect) + c.B;
    case 'exp_1_0'
        FitCurve = exp(-c.n*Dvect);
    case 'decay'
        FitCurve = c.A/Dvect + c.B;
        Residuals = Rvect' - FitCurve;
    case 'exp0'
        FitCurve = c.A.*exp(-c.n*Dvect);
    case 'exp1'
        FitCurve = exp(-c.n*Dvect) + c.B;
    otherwise
        Y = discretize(distExpVect(:,1),xThresholds);
        Residuals = zeros(length(Y),1);
        for val=1:length(Y)
            if ~isnan(distExpVect(val,2))
                Residuals(val) = distExpVect(val,2) - yMeans(Y(val));
            else
                Residuals(val) = NaN;
            end
        end
end
end