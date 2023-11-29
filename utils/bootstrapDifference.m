function [sigdeltaPD, sigdeltaMD] = bootstrapDifference(tunings, refDay, n, alpha, seed)

% check if the confidence intervals of the bootstrap difference distribution
% between two distributions of tuning contains zero
% If the 95% confidence interval for the difference distribution does not contain 0, 
% then we reject the null hypothesis at the 5% significance level. 
% adopted method from  I. H. Stevenson et al., 
% "Statistical assessment of the stability of neural movement representations,” 
% J. Neurophysiol., vol. 106, pp. 764–774, 2011, [Online]. 

if nargin < 5; seed = 42; end
if nargin < 4; alpha = 0.05; end
if nargin < 3; n = 1000; end
rng(seed)

[nfeats, nday] = size(tunings.pd);
sigdeltaPD = NaN(nfeats, nday);
sigdeltaMD = NaN(nfeats, nday);

% bootstrap to get confidence intervals
for f = 1:nfeats
    for day = refDay(f):nday
        if (tunings.pval(f,day) >= alpha)
            continue
        else
            % PD
            perm1 = randi(n,1,n); perm2 = randi(n,1,n);
            pd_1 = tunings.bootstrap.pd(f,perm1,refDay(f))';
            pd_2 = tunings.bootstrap.pd(f,perm2,day)';
            deltaPDdist = angdiffdeg(pd_1,pd_2);
            pd_mu = circ_mean(deltaPDdist*pi/180)*180/pi;
            pd_ci = prctile(deltaPDdist,[100*alpha/2,100*(1-alpha/2)]);
            sigdeltaPD(f,day) = all(pd_ci>0) || all(pd_ci<0);
    
            % MD
            md_1 = tunings.bootstrap.md(f,perm1,refDay(f))';
            md_2 = tunings.bootstrap.md(f,perm2,day)';
            deltaMDdist = (md_2 - md_1);
            md_mu = circ_mean(deltaMDdist*pi/180)*180/pi;
            md_ci = prctile(deltaMDdist,[100*alpha/2,100*(1-alpha/2)]);
            sigdeltaMD(f,day) = all(md_ci>0) || all(md_ci<0);
        end
    end
end

% make example histogram
figure()
subplot(1,2,1)
histogram(deltaPDdist,(-180:5:180))
histogram(pd_1,(-180:5:180))
histogram(pd_2,(-180:5:180))
xline(pd_ci);xline(pd_mu)
subplot(1,2,2)
bins = linspace(-0.1,0.1,101);
histogram(deltaMDdist,bins)
histogram(md_1,bins)
histogram(md_2,bins)
xline(md_ci);xline(md_mu)

end