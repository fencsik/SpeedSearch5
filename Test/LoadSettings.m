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
end

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
    end
end

function fid = FindAndOpenSettingsFile ()
    global par;
    found = 0;
    for i = 1:numel(par.settingsFileName)
        fileName = par.settingsFileName{i};
        if exist(fileName, 'file')
           found = 1;
           break;
        end
    end
    if found
        fid = fopen(fileName, 'r');
        fprintf('Opened settings file %s\n', fileName);
    else
        fid = -1;
        fprintf('Could not open settings file\n');
    end
end

function settingsStruct = ProcessSettingsFile (fid)
    settingsStruct = struct();
    key = 1;
    while (1)
        [key, value] = GetNextTokenFromSettingsFile(fid);
        if (isempty(key))
            break;
        else
            settingsStruct = setfield(settingsStruct, key, value);
        end
    end
end

function [key, val] = GetNextTokenFromSettingsFile (fid)
    key = [];
    val = [];
    while (1)
        str = fgetl(fid);
        if (str == -1)
            break;
        end
        str = strtrim(RemoveComments(str));
        if (~isempty(str))
            tokens = strtrim(strsplit(str, '='));
            if (numel(tokens) >= 1)
                key = ProcessKey(tokens{1});
            end
            if (numel(tokens) >= 2)
                val = ProcessValue(tokens{2});
            end
            break;
        end
    end
end

function key = ProcessKey(keyToken)
    key = StripWhitespaceFromString(keyToken);
end

function value = ProcessValue(valueToken)
    value = SplitStringIntoVector(valueToken);
    value = ConvertToNumbers(value);
    value = ConvertToVector(value);
end

function stringOut = RemoveComments(stringIn)
    % remove comment character and everything after it
    stringOut = regexprep(stringIn, '#.*', '');
end

function stringOut = StripWhitespaceFromString (stringIn)
    % remove all whitespace from a string
    stringOut = regexprep(stringIn, '\\s', '');
end

function vector = SplitStringIntoVector(stringIn)
%%% if stringIn has separate components, split them
    %% replace array indicators/separators with whitespace
    stringIn = regexprep(stringIn, '\\[|\\]|,', ' ');
    stringIn = strtrim(stringIn);
    vector = strsplit(stringIn);
    if (numel(vector) == 1)
        vector = vector{1};
    end
end

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
        end
    else
        for (i = 1:numel(vectorIn))
            n = str2double(vectorIn{i});
            if (isnan(n)) % conversion to number failed
                vectorOut{i} = vectorIn{i};
            else
                vectorOut{i} = n;
            end
        end
    end
end

function vectorOut = ConvertToVector(vectorIn)
%%% Converts numeric cell arrays to vectors.  Leaves cell arrays with
%%% strings alone.
    if (iscell(vectorIn) && IsCellArrayNumeric(vectorIn))
        vectorOut = ConvertNumericCellArrayToVector(vectorIn);
    else
        vectorOut = vectorIn;
    end
end

function val = IsCellArrayNumeric(vector)
    if (iscell(vector))
        val = 1;
        for (i = 1:numel(vector))
            if (ischar(vector{i}))
                val = 0;
                break;
            end
        end
    else
        val = 0;
    end
end

function out = ConvertNumericCellArrayToVector(in)
    n = numel(in)
    out = zeros(1, n)
    for (i = 1:n)
        out(i) = in{i};
    end
end

function MakeSettingsGlobal (settingsArray)
    global par
    fn = fieldnames(settingsArray);
    for (i = 1:numel(fn))
        par.(fn{i}) = settingsArray.(fn{i});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Local Variables:
%%% mode:Octave
%%% End:
