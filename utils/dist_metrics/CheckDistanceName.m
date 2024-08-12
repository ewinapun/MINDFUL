function [dfnNew, dfnStr] = CheckDistanceName(dfn)

% Ensure distance names are correctly spelled
% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

    alldfn = {'Jdiv','meanM','WcovM','KLdiv','Bhatt','Wass'};
    alldfnstr = {'Jeffrey''s divergence (symm KL)',...
               'Mean difference',...
               'Covariance (Wasserstein)', ...
               'KL divergence', ...
               'Bhattacharyya', ...
               'Wasserstein'};
    if isempty(dfn) 
        dfnNew = alldfn;
        dfnStr = alldfnstr;
        sp1 = 2;
        sp2 = 3;
    else
        if ~iscell(dfn);dfn = cellstr(dfn);end 
        dfnNew = [];
        dfnStr = [];
        dfn = unique(lower(dfn),'stable');
        D = length(dfn);
        for d = 1:D
            loc = strcmpi(dfn{d},alldfn);
            dfnNew = [dfnNew alldfn(loc)];
            dfnStr = [dfnStr alldfnstr(loc)];
        end
    end
end