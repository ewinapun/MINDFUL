function out = CallAffineWarp(features, params)
% features should be nTrials x nTime x nDim
% See dfParams for defaults and bgAffineWarp.py for all params.
% 
% See also PlotWarpParams ComputeWarpTrialVal

saveFile = 'E:\Data\affine_warp\dpcaWarpAlignIn.mat';
pyOutputFileLoc = 'E:\Data\affine_warp\dpcaWarpAlignOut.mat';
batchFileLoc = 'C:\Users\Tommy\Code\Repos\analysis\UserAnalysis\hosman\Python\tutils\affineDayToStartup.bat ';
if nargin == 0
    out = load( pyOutputFileLoc );
    return
end

if ndims(features) ~= 3
    error('Expected features to have 3 dim, found %d\n', ndims(features))
end
if nargin < 2
    params = [];
end
% Default params are in bgAffineWarp.py > BGAffineWarp > get_default_params
dfParams.warp_reg_scale = 0.001; % Trial warp
dfParams.smoothness_reg_scale = 0.3; % Feature Template smoothness
dfParams.l2_reg_scale = 1e-3; % Feature Template 
dfParams.fit_iterations = 15;
dfParams.warp_iterations = 50;
dfParams.saveLocation = pyOutputFileLoc; % Where the python program saves its output
dfParams.save_plot_location = fullfile(fileparts(pyOutputFileLoc),'dPCAWarpPlots');
dfParams.model_knots = [-1 0]; %[-1 0 1]; % [-1 0 1 2]; % [-1,0];
params = MergeBintoA(dfParams,params);


%% Call affine warp


warning off
try; delete(params.saveLocation); end
warning on


fprintf('Saving...\n%s\n',saveFile)
save( saveFile, 'features', 'params' );
fprintf('Done!\n')
tic
system([batchFileLoc saveFile])
toc

saveFeatWhos = whos('features');
fprintf('Waiting for %s to be saved...\n', params.saveLocation);

out = WaitForFile(params, saveFeatWhos);



end


function out = WaitForFile(params, saveFeatWhos)


aa = dir(params.saveLocation);
while isempty(aa) || aa.bytes < saveFeatWhos.bytes
    pause(0.5)
    aa = dir(params.saveLocation);
end


try
    % Load warped_data
    out = load( params.saveLocation );
    fprintf('Loading complete!\n')
catch err
    beep
    UnrollError(err);
    fprintf(2, '\n\nKeyboard at error\n')
    keyboard
end
    

end