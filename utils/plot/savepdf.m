function savepdf(gcf, filename, opt)
% helper function to save figures into pdf
% inputs:
%     gcf       - Get handle to current figure
%     filename  - filename of the figure to be saved
%     opt       - 'fig': optional to also save as fig
%               - 'png': optional to also save as png 

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

    if nargin < 3
        opt = '';
    end
    
    % check if filename already exist to prevent overwriting
    currentFolder = cd;
    while (exist([fullfile(currentFolder, filename) '.pdf'], 'file') == 2)
        rename = input('Warning: file already exist, do you want to rename?[y/n]','s');
        if rename == 'y'
            filename = input('New filename?','s');
        else
            break
        end
    end
    % save as pdf
    set(gcf,'Units','inches');
    set(gcf,'Color','white');
    set(gcf, 'InvertHardcopy', 'off');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print(gcf, filename, '-dpdf', '-loose')
    disp(['Saving image as ',filename])
    if strcmp(opt, 'fig')
        saveas(gcf, filename,'fig')
    end
    if strcmp(opt, 'png')
        saveas(gcf, filename,'png')
    end
end