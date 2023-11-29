function sTo = ConcatStruct(sTo, sFrom, dim)
% sTo = ConcatStruct(sTo, sFrom, dim)
% sFrom -> sTo 
% along dim if passed otherwise along the first non-singleton.
 

if nargin < 3
    dim = [];
end

try
    if isempty(sTo)
        sTo = sFrom;
    else
        fn = fieldnames(sTo);
        for jj = 1:length(fn)
            
            if isstruct(sTo.(fn{jj}))
                % Recurse
                sTo.(fn{jj}) = ConcatStruct(sTo.(fn{jj}), sFrom.(fn{jj}), dim);
            else
                
                if isempty(dim)
                    catDim = GetDim(sTo.(fn{jj}),sFrom.(fn{jj}));
                else
                    catDim = dim;
                end
                
                
                % Concat
                sTo.(fn{jj}) = cat(catDim, sTo.(fn{jj}), sFrom.(fn{jj}));
            end
            
        end
    end
catch err
    UnrollError(err)
    keyboard
    rethrow(err)
end
end

function dim = GetDim(to,from)
    % First non-singleton dimension, but ignore empty
    % that is, > 1
    catDimTo = find( size(to) > 1,1);
    if isempty(catDimTo)
        catDimTo = nan;
    end
    
    catDimFrom = find( size(from) > 1,1);
    if isempty(catDimFrom)
        catDimFrom = nan;
    end
    
    dim = nanmin(catDimTo,catDimFrom);
    
    if isnan(dim)
        dim = 1;
    end
end