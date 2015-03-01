function LoadSettings ()
    # Tests code to load settings from a settings file.  Intended to be put
    # directly into Psychtoolbox scripts.
    global par = struct();
    par.settingsFileName = {"Settings.txt", "Settings", ...
                            "settings.txt", "settings"};
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

####################################################

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
        if exist(fileName, "file")
           found = 1;
           break;
        endif
    endfor
    if found
        fid = fopen(fileName, "r");
        printf("Opened settings file %s\n", fileName);
    else
        fid = -1;
        printf("Could not open settings file\n");
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
            tokens = strtrim(strsplit(str, "="));
            key = tokens{1};
            val = tokens{2};
            break;
        endif
    endwhile
endfunction

function stringOut = RemoveComments(stringIn)
    # remove comment character and everything after it
    stringOut = regexprep(stringIn, "#.*", "");
endfunction

function stringOut = StripWhitespaceFromString (stringIn)
    # remove all whitespace from a string
    stringOut = regexprep(stringIn, "\\s", "");
endfunction

function MakeSettingsGlobal (settingsArray)
    global par
    fn = fieldnames(settingsArray);
    for (i = 1:numel(fn))
        par.(fn{i}) = settingsArray.(fn{i});
    endfor
endfunction

###############################################


### Local Variables:
### mode:Octave
### End:
