function errorStrOut = UnrollError( sError, useHyperlinks, printErr )
% This function unrolls and prints an error structure from a try catch
% statement.
%
% Example
%
% try
%    code...
% catch sError 
%    % Print out error message
%    errorStr = UnrollError( sError );
% 
% end
%
% Written by Tommy Hosman


% useHyperlinks in the output (this is so you can click on the error in the
% matlab console)
if ~exist('useHyperlinks','var') || isempty(useHyperlinks)
    useHyperlinks = 1;
end
if ~exist('printErr','var') || isempty(printErr)
    printErr = 1;
end

    
    errorStr = sprintf('\nERROR - %s\n\n', sError.message);
    
    for stackInd = 1:length(sError.stack)
        [filePath,fileName] = fileparts(sError.stack(stackInd).file);
        fullFunction = [filePath filesep fileName];
        
        lineNum = sError.stack(stackInd).line;
        
        if useHyperlinks
            lineNuberStr = sprintf('<a href="matlab: opentoline(''%s.m'',%d,0)">line %d</a>', fullFunction, lineNum, lineNum);
        else
            lineNuberStr = sprintf('line %d', lineNum);
        end
        
        errorStr = [errorStr ...
                    sprintf('file %s - function %s (%s)\n',...
                            fileName, sError.stack(stackInd).name, lineNuberStr)];
                        
    end
    
    if printErr
        % Print as error
        fprintf(2, '%s\n', errorStr );
    end
    
    
    if nargout > 0
        errorStrOut = errorStr;
    end
    
    
end