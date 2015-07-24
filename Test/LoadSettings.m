function LoadSettings ()
    % Tests code to load settings from a settings file.  Intended to be put
    % directly into Psychtoolbox scripts.
    global par = struct();
    par.settingsFileName = {'Settings.txt', 'Settings', ...
                            'settings.txt', 'settings'};
    disp(par);
    LoadSettingsFile();
    disp(par);
    clear -global
endfunction

function success = CloseFile(fid)
    x = fclose(fid);
    if (x == 0)
        success = 1;
    else
        success = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LoadSettingsFile ()
    fid = FindAndOpenSettingsFile();
    if fid != -1
        settingsArray = ProcessSettingsFile(fid);
        MakeSettingsGlobal(settingsArray);
        CloseFile(fid);
    endif
endfunction

function fid = FindAndOpenSettingsFile ()
    global par;
    found = 0;
    for i = 1:numel(par.settingsFileName)
        fileName = par.settingsFileName{i};
        if exist(fileName, 'file')
           found = 1;
           break;
        endif
    endfor
    if found
        fid = fopen(fileName, 'r');
        printf('Opened settings file %s\n', fileName);
    else
        fid = -1;
        printf('Could not open settings file\n');
    endif
endfunction

function settingsStruct = ProcessSettingsFile (fid)
    settingsStruct = struct();
    key = 1;
    while (1)
        [key, value] = GetNextTokenFromSettingsFile(fid);
        if (isempty(key))
            break;
        else
            settingsStruct = setfield(settingsStruct, key, value);
        endif
    endwhile
endfunction

function [key, val] = GetNextTokenFromSettingsFile (fid)
    key = [];
    val = [];
    while (1)
        str = fgetl(fid);
        if (str == -1)
            break;
        endif
        str = strtrim(RemoveComments(str));
        if (~isempty(str))
            tokens = strtrim(strsplit(str, '='));
            if (numel(tokens) >= 1)
                key = ProcessKey(tokens{1});
            endif
            if (numel(tokens) >= 2)
                val = ProcessValue(tokens{2});
            endif
            break;
        endif
    endwhile
endfunction

function key = ProcessKey(keyToken)
    key = StripWhitespaceFromString(keyToken);
endfunction

function value = ProcessValue(valueToken)
    value = SplitStringIntoVector(valueToken);
    value = ConvertToNumbers(value);
    value = ConvertToVector(value);
endfunction

function stringOut = RemoveComments(stringIn)
    % remove comment character and everything after it
    stringOut = regexprep(stringIn, '#.*', '');
endfunction

function stringOut = StripWhitespaceFromString (stringIn)
    % remove all whitespace from a string
    stringOut = regexprep(stringIn, '\\s', '');
endfunction

function vector = SplitStringIntoVector(stringIn)
%%% if stringIn has separate components, split them
    %% replace array indicators/separators with whitespace
    stringIn = regexprep(stringIn, '\\[|\\]|,', ' ');
    stringIn = strtrim(stringIn);
    vector = strsplit(stringIn);
    if (numel(vector) == 1)
        vector = vector{1};
    endif
endfunction

function vectorOut = ConvertToNumbers(vectorIn)
%%% Converts number strings to floats.  Leaves character strings alone
    if (isempty(vectorIn))
        vectorOut = vectorIn;
    elseif (!iscell(vectorIn))
        n = str2double(vectorIn);
        if (isnan(n))
            vectorOut = vectorIn;
        else
            vectorOut = n;
        endif
    else
        for (i = 1:numel(vectorIn))
            n = str2double(vectorIn{i});
            if (isnan(n)) % conversion to number failed
                vectorOut{i} = vectorIn{i};
            else
                vectorOut{i} = n;
            endif
        endfor
    endif
endfunction

function vectorOut = ConvertToVector(vectorIn)
%%% Converts numeric cell arrays to vectors.  Leaves cell arrays with
%%% strings alone.
    if (iscell(vectorIn) && IsCellArrayNumeric(vectorIn))
        vectorOut = ConvertNumericCellArrayToVector(vectorIn);
    else
        vectorOut = vectorIn;
    endif
endfunction

function val = IsCellArrayNumeric(vector)
    if (iscell(vector))
        val = 1;
        for (i = 1:numel(vector))
            if (ischar(vector{i}))
                val = 0;
                break;
            endif
        endfor
    else
        val = 0;
    endif
endfunction

function out = ConvertNumericCellArrayToVector(in)
    n = numel(in)
    out = zeros(1, n)
    for (i = 1:n)
        out(i) = in{i};
    endfor
endfunction

function MakeSettingsGlobal (settingsArray)
    global par
    fn = fieldnames(settingsArray);
    for (i = 1:numel(fn))
        par.(fn{i}) = settingsArray.(fn{i});
    endfor
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Local Variables:
%%% mode:Octave
%%% End:
