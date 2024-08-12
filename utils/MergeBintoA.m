function a = MergeBintoA(a, b, options)
% a = MergeBintoA(a, b, options)
% 
% Recursively puts all fields from b to a.
% 
% 
% • Values that are unique to a will remain in a. 
% • Values that are unique to b will be in a. 
% • If a and b have the same substructure, then b will overwrite a's value.
% • If a and b have the same substructure name of different types and
%   strictTypes is set true, the function will error. See CheckTypes()
% 
% This can be used to set default parameters (as struct 'a') and
% user-defined parameters (as struct 'b') to allow a user to overwrite the
% default params.
% 
% Note! The substructure fieldname order may change after using this
% structure.
% 
% 
% Inputs:
% 
%   a 
%       structure with the default fields, to be replaced or merged
% 
%   b
%       contains fields we want to copy to a
% 
%   options
% 
%       .copyDirection (default: 'auto')
%           A string input specify the direction to copy the merged
%           structures copyDirection can be: 'b2a', 'a2b', 'auto' default:
%           'auto'
% 
%           Fields in b will always override a, but you may choose to copy
%           from a into b or from b into a depending on struct size.
% 
%           Auto performs 'whos' on both structures and copies fields of
%           the smaller structure to the larger one, again, with fields in
%           struct b overriding fields of struct a.
%       
%       .strictTypes (default: false)
%           Check that a and b have the same object type before copying
%           (e.g. are they both 'double'?).
% 
%       .noMergeIfEmpty (default: false)
%           Some instances, a field exists in structure b, but if it is
%           empty, we may want to have the value in struct a remain.
%
% Outputs:
%   a - structure that has its values replaced and merged with struct b
%
% 
% Example: (note subfield a.ee gets overwritten)
% 
%     a.a    = 1;
%     a.b    = 2;
%     a.c    = 3;
%     a.ee.e = 123;
% 
%     b.b    = 1;
%     b.d    = 10;
%     b.ee   = 12;
% 
%     a = MergeBintoA(a,b)
% 
%     a = 
% 
%          a: 1
%          b: 1
%          c: 3
%          d: 10
%         ee: 12
% 
% 
% 
% Example 2: With strictTypes = true
% 
%     a.a    = 1;
%     a.b.b  = 2;
%     a.c    = 3;
%     a.ee.e = 123;
% 
%     b.b    = 1;
%     b.d    = 10;
%     b.ee   = 12;
% 
%     options.strictTypes = true;
%     a = MergeBintoA(a,b, options)
% 
%     Error using MergeBintoA>CheckTypes (line 433)
%     Types do not match between values for field ee. Input is type 'double' and output is of type 'struct'
% 
%     Error in MergeBintoA (line 112)
%                 CheckTypes( b, a, inFN{i} )
% 
% 
% 
% 
% History:
%   2017   Copyright Tommy Hosman, Brown University. All Rights Reserved
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Handle Inputs
%--------------------------------------------------------------------------
if nargin < 3
    options = [];
end
if ~isempty(options) && ~isstruct(options)
    error('options must be a structure. You passed in a variable of type %s', class(options))
end

if ~isfield(options, 'copyDirection')
    options.copyDirection = 'auto';
end
if ~isfield(options, 'strictTypes')
    options.strictTypes = 0;
end
if ~isfield(options, 'printOn')
    options.printOn = 0;
end
if ~isfield(options, 'noMergeIfEmpty')
    options.noMergeIfEmpty = 0;
end

% Do nothing if b is not a struct
if ~isstruct(b)
    return;
end

% If a is not a struct, set b to a
if ~isstruct(a)
    a = b;
    return;
end



%--------------------------------------------------------------------------
% Auto detect copy direction
%--------------------------------------------------------------------------
if strcmpi(options.copyDirection,'auto')
    aInfo = whos('a');
    bInfo = whos('b');
    if aInfo.bytes > bInfo.bytes
        % Larger a struct, copy smaller b into a
        options.copyDirection = 'b2a';
    else
        % Larger b struct, copy smaller a into fields of b that do not exist
        options.copyDirection = 'a2b';
    end
end




%--------------------------------------------------------------------------
%  Get variable names for printing, messages and errors
%--------------------------------------------------------------------------
rootName = {inputname(1), inputname(2)};

if options.printOn
    fprintf('\nMerging %s -> %s\n\n', rootName{2}, rootName{1} );
end



%--------------------------------------------------------------------------
%  Merge structures
%--------------------------------------------------------------------------
switch options.copyDirection
    case 'b2a'
        mergeAll = 1;
        a = MergeStructsRecursively(b, a, mergeAll, options.strictTypes, rootName, options.noMergeIfEmpty, options.printOn);
    case 'a2b'
        mergeAll = 0;
        a = MergeStructsRecursively(a, b, mergeAll, options.strictTypes, rootName, options.noMergeIfEmpty, options.printOn);
        
        
    otherwise
        error('You have specified an invalid option ''%s''', options.copyDirection )
end



if options.printOn
    fprintf('Merge complete\n\n');
end

end









%% Recursive merge function

function s2 = MergeStructsRecursively(s1, s2, mergeAll, strictTypes, rootName, noMergeIfEmpty, printOn)
% function s2 = MergeStructsRecursively(s1, s2, mergeAll, strictTypes, rootName)
% 
% Merges structure s1 -> s2. 
% 
% If mergeAll is set true, all subfields of struct 's1' are copied to 's2'.
% If mergeAll is set false, only fields that are unique to struct s1 are
% copied to s2. See below for more details and examples.
% 
% If strictTypes is set true, this function will error if copied fields are
% of different types.
% 
% 
% 
% Inputs:
%   s1
%       Structure whose subfields are recursively copied to stucture s2
%   s2
%       Structure where structure s1's subfields are copied to.
% 
%   mergeAll
%       Flag that indicates what fields of structure a are copied to
%       structure b. If true, all fields are copied from a to b. If false,
%       only fields that are not common between a and b are copied. See
%       below for more details.
% 
%   strictTypes
%       Flag that if set true will only copy fields that are of the same
%       class (e.g. double, or uint8). If two fields are of different
%       classes, and strictTypes is true, then this function will raise an
%       error.
% 
%       If false, no type checking is performed.
% 
%   rootName
%       A cell array consisting of the names of each structure passed in.
% 
%   noMergeIfEmpty
%       Flag that if set true, we will overwrite b, if its subfield is
%       empty. If false, a subfield that is empty will not be overwritten.
% 
%   printOn
%       Flag that if set true will print the fields that are overwritten.
% 
% 
% 
% Output
% 
%   s2
%       Structure that has all the merged fields.
% 
% 
% -More details on input flag "mergeAll"
% 
% If mergeAll is true:
% 
%   s2 = s2 <- s1
% 
%   All fields from struct 's1' will be pushed to struct 's2'. But any
%   unique fields of struct b will remain in s2.
% 
%   Put another way, sub-fields of struct s2 that are common with struct
%   s1, will be overwritten by the data in struct s1.
%   
% 
% If mergeAll is false:
%   
%   s2 = s1 <- setdiff(s1,s2) 
% 
%   Only fields that are not in struct 's2' will be copied from struct 's1'
%   to struct 's2'. All fields that are common between s2 and s1 will keep
%   the values from struct s2.
% 
%   That is, sub-fields that are common with struct s2 and s1 will NOT be
%   overwritten by the data in struct s1.
% 
% 
% 
% 
% Example 1 - mergeAll is true
%     s1.a    = 1;
%     s1.b    = 2;
%     s1.c    = 3;
%     s1.ee.e = 123;
% 
%     s2.b    = 1;
%     s2.d    = 10;
%     s2.ee   = 12;
% 
%     strictTypes = false;
%     mergeAll = 1;
%     s2 = MergeStructsRecursively(s1, s2, mergeAll, strictTypes)
% 
%     % output:
%     s2 = 
% 
%       struct with fields:
% 
%          b: 2
%          d: 10
%         ee: [1×1 struct]
%          a: 1
%          c: 3
% 
% 
% Example 2 - mergeAll is false
%     s1.a    = 1;
%     s1.b    = 2;
%     s1.c    = 3;
%     s1.ee.e = 123;
% 
%     s2.b    = 1;
%     s2.d    = 10;
%     s2.ee   = 12;
% 
%     strictTypes = false;
%     mergeAll = 0;
%     s2 = MergeStructsRecursively(s1, s2, mergeAll, strictTypes)
% 
%     % output:
%     s2 = 
% 
%       struct with fields:
% 
%          b: 1
%          d: 10
%         ee: 12
%          a: 1
%          c: 3
% 
% 
% 

%% Code



% Get field names for structure 's1'
s1_FN = fieldnames(s1);


%--------------------------------------------------------------------------
% Loop through all fields and attempt to set them to the output struct (s2)
%--------------------------------------------------------------------------
for ii = 1:length( s1 )
    
for i = 1:length(s1_FN)
    
    
    
    %----------------------------------------------------------------------
    % Append the subfield name to the root
    %----------------------------------------------------------------------
    fullFieldName{1} = [rootName{1} '.' s1_FN{i}];
    fullFieldName{2} = [rootName{2} '.' s1_FN{i}];
    
    
    if ii > length( s2 )
        isFieldInS2 = false;
        s2 = MergeData(s1, s2, s1_FN, i, ii, mergeAll, noMergeIfEmpty, isFieldInS2, fullFieldName, printOn);
        continue;
    end
    
    
    %----------------------------------------------------------------------
    % Is this field a substructure in both 's1' and 's2'
    %----------------------------------------------------------------------
    isFieldInS2 = isfield( s2(ii), s1_FN{i} );
    if isstruct( s1(ii).(s1_FN{i}) ) && isFieldInS2 && isstruct( s2(ii).(s1_FN{i}) )
        
        
        %--------------------------------------------------------------------------
        % Recursively call ourselves to parse the substructure
        %--------------------------------------------------------------------------
        s2(ii).(s1_FN{i}) = MergeStructsRecursively( s1(ii).(s1_FN{i}), s2(ii).(s1_FN{i}),  mergeAll, strictTypes, rootName, noMergeIfEmpty, printOn );
        
        
    
        
    else
        
        %------------------------------------------------------------------
        % Here the subfield is not a struct for s1. Merge.
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % If strictTypes is true, confirm that the two fields are from
        % the same class.
        %------------------------------------------------------------------
        if strictTypes
            CheckTypes( s1(ii).(s1_FN{i}), s2.(s1_FN{i}), fullFieldName );            
        end
        
        
        s2 = MergeData(s1, s2, s1_FN, i, ii, mergeAll, noMergeIfEmpty, isFieldInS2, fullFieldName, printOn);
       
        
        
        
        
        
        
    end
    
end
end

% Remove elements of s2 that were not in s1
if length(s2) > length(s1)
    s2(end+1:end) = [];
end

end

function s2 = MergeData(s1, s2, s1_FN, i, ii, mergeAll, noMergeIfEmpty, isFieldInS2, fullFieldName, printOn)
 %------------------------------------------------------------------
% If noMergeIfEmpty is true, we will want to push s1 into s2 if s1
% is not empty.
%------------------------------------------------------------------
if noMergeIfEmpty && isFieldInS2
    mergeIfNotEmpty = ~isempty( s1(ii).(s1_FN{i}) );
else
    mergeIfNotEmpty = true;
end
%------------------------------------------------------------------
% Merge if mergeAll was requested, or if the field does not exist.
%------------------------------------------------------------------        
% That is, if mergeAll was requested, always push from s1 to s2
% If not, only push from s1 to s2 if the field does not exist in s2.       

if ( mergeAll || ~isFieldInS2 ) && mergeIfNotEmpty


    if printOn
        fprintf('Updating \t %s\n', fullFieldName{2});
    end


    %--------------------------------------------------------------
    % Merge here
    %--------------------------------------------------------------

    s2(ii).(s1_FN{i}) = s1(ii).(s1_FN{i});



end    

end

%--------------------------------------------------------------------------
% CheckTypes
%--------------------------------------------------------------------------
function CheckTypes( s1, s2, abFullfieldnames )
% Make sure that both s1 and s2 are of the same type (e.g. double or uint8)

 
    % Do the two objects have the same type?
    if ~strcmpi(class(s1), class(s2))
        
        error(['Types do not match between values for subfield\n', ...
               '%s \tis type ''%s''\n', ...
               '%s \tis type ''%s'''], ...           
               abFullfieldnames{1}, class(s1), ...
               abFullfieldnames{2}, class(s2))
    end

end
