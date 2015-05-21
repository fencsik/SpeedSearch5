function varargout = Staircaser (command, varargin)
    global _staircaserPar;
    if (!IsStaircaserInitialized())
        InitializeStaircaser();
    endif
    _staircaserPar.nOutputArgs = nargout;

    _staircaserPar.printHelp = 0;
    if command(end) == "?"
        _staircaserPar.printHelp = 1;
        command = command(1:end-1);
    endif

    switch command
        case "Create"
            varargout = StaircaserCreate(varargin);
        case "Delete"
            varargout = StaircaserDelete(varargin);
        case "StartTrial"
            varargout = StaircaserStartTrial(varargin);
        case "EndTrial"
            varargout = StaircaserEndTrial(varargin);
        case "FinalValue"
            varargout = StaircaserFinalValue(varargin);
        case "GetReversals"
            varargout = StaircaserGetReversals(varargin);
        case "Plot"
            varargout = StaircaserPlot(varargin);
        case "Progress"
            varargout = StaircaserProgress(varargin);
        case "List"
            varargout = StaircaserList(varargin);
###     case "PrintTracks"
###         varargout = StaircasePrintTracks(varargin);
###     case "GetValue"
###         varargout = StaircaserGetValue(varargin);
###     case "GetTrack"
###         varargout = StaircaserGetTrack(varargin);
###     case "IsDone"
###         varargout = StaircaserIsDone(varargin);
        otherwise
            varargout = [];
            error("Staircaser: command %s not recognized", command);
    endswitch
endfunction

function InitializeStaircaser ()
    global _staircaserPar;
    _staircaserPar = struct();
    _staircaserPar.debug = 2;
    _staircaserPar.nTrialsBase = 4;
    _staircaserPar.staircase = struct();
    _staircaserPar.idList = [];
    if (_staircaserPar.debug >= 1)
        printf("%s: Initialized staircaser\n", mfilename);
    endif
endfunction

function tf = IsStaircaserInitialized()
    global _staircaserPar;
    tf = 0;
    if (exist("_staircaserPar", "var") && isstruct(_staircaserPar) &&
        isfield(_staircaserPar, "staircase") && isstruct(_staircaserPar.staircase))
        tf = 1;
    endif
endfunction

function id = FindUnusedStaircaseID ()
    global _staircaserPar;
    idList = _staircaserPar.idList;
    if (isempty(idList))
        id = 1;
        return;
    endif
    ## find the first unused number
    for (i = 1:max(idList))
        if (!any(i == idList))
            id = i;
            return;
        endif
    endfor
    ## if we've gotten here, then all numbers up to max(idList) are in use
    id = max(idList) + 1;
endfunction

function CheckNumberOfInputArguments(nArgs, minArgs, maxArgs, helpText=[])
    if (nArgs < minArgs)
        error("%s\nNot enough input arguments", helpText);
    elseif (nArgs > maxArgs)
        error("%s\nToo many input arguments", helpText);
    endif
endfunction

function CheckNumberOfOutputArguments(nArgs, minArgs, maxArgs, helpText=[])
    if (nArgs < minArgs)
        error("%s\nNot enough output arguments", helpText);
    elseif (nArgs > maxArgs)
        error("%s\nToo many output arguments", helpText);
    endif
endfunction

function helpRequested = CheckForHelpRequest (helpText)
    global _staircaserPar;
    helpRequested = 0;
    if (_staircaserPar.printHelp)
        printf('%s', helpText);
        helpRequested = 1;
    endif
endfunction

###########################################################################
### Staircaser("Create")

function argout = StaircaserCreate(argin)
    global _staircaserPar;
    helpText = sprintf(["\nid = Staircaser(\"Create\", type, nReversals, initial,\n", ...
                        "steps, [nReversalsDropped], [nTracks], [range]);\n\n"]);
    if (CheckForHelpRequest(helpText))
        argout = {[]};
        return;
    endif

    ## process input arguments
    nargs = numel(argin);
    CheckNumberOfInputArguments(nargs, 4, 7, helpText);
    type = argin{1}(1);
    if (type != 1)
        error("staircase types other than 1 are not yet supported");
    endif
    nReversals = argin{2}(1);
    initialValue = argin{3}(1);
    steps = argin{4};
    if nargs < 5 || isempty(argin{5})
        nReversalsDropped = 0;
    else
        nReversalsDropped = argin{5}(1);
    endif
    if nargs < 6 || isempty(argin{6})
        nTracks = 1;
    else
        nTracks = argin{6}(1);
        nTracks = nTracks(1);
    endif
    if nargs < 7 || numel(argin{7}) < 2
        range = [];
    else
        range = [min(argin{7}), max(argin{7})];
    endif
    ## process output arguments
    CheckNumberOfOutputArguments(_staircaserPar.nOutputArgs, 0, 1, helpText);

    ## error handling
    if !any(sign(steps) > 0)
        error("no positive steps");
    elseif ~any(sign(steps) < 0)
        error("no negative steps");
    elseif nReversals == nReversalsDropped
        error("nReversalsDropped = nReversals, so no reversals would be used");
    endif

    staircase = struct();
    staircase.type = type;
    staircase.nReversals = nReversals;
    staircase.steps = steps;
    staircase.nReversalsDropped = nReversalsDropped;
    staircase.nTracks = nTracks;
    staircase.nTracksRemaining = nTracks;
    staircase.range = range;
    staircase.inTrial = 0;
    staircase.currentTrack = 0;
    staircase.lastTrack = 0;
    staircase.isDone = 0;
    staircase.progress = 0;
    staircase.progressTotal = nTracks * nReversals;
    staircase.finalValue = [];
    staircase.reversals = [];
    staircase.values = [];
    staircase.responses = [];
    scLabels = cell(nTracks, 1);
    for i = 1:nTracks
        scLabels{i} = i;
    endfor
    staircase.tracks = struct("label", scLabels,
                              "value", initialValue,
                              "counter", 0,
                              "lastStepDir", [],
                              "reversals", nan(nReversals, 1),
                              "values", nan(_staircaserPar.nTrialsBase, 1),
                              "responses", nan(_staircaserPar.nTrialsBase, 1),
                              "nTrials", 0);
    ## get next available id
    id = FindUnusedStaircaseID();
    _staircaserPar.idList = [_staircaserPar.idList; id];
    _staircaserPar.staircase(id) = staircase;
    if _staircaserPar.debug >= 1
        fprintf(["%s: created staircase %d with %d reversals, %d ", ...
                 "reversals dropped, %d tracks\n"], mfilename, id, ...
                nReversals, nReversalsDropped, nTracks);
        fprintf("%s: and range of \n", mfilename);
        disp(range);
    endif
    if _staircaserPar.debug >= 2
        for i = 1:nTracks
            disp(staircase.tracks(i));
        endfor
    endif
    argout = {id};
endfunction

###########################################################################
### Staircaser("Delete")

function argout = StaircaserDelete(argin)
    global _staircaserPar;
    helpText = sprintf("\nsuccess = Staircaser(\"Delete\", id);\n\n");
    if (CheckForHelpRequest(helpText))
        argout = {[]};
        return;
    endif

    ## process input arguments
    nargs = numel(argin);
    CheckNumberOfInputArguments(nargs, 1, 1, helpText);
    ## process output arguments
    CheckNumberOfOutputArguments(_staircaserPar.nOutputArgs, 0, 1, helpText);

    id = argin{1};
    idList = _staircaserPar.idList;
    success = zeros(1, numel(id));
    for i = 1:numel(id)
        if (!isnumeric(id(i)))
            error("invalid staircase id given");
        elseif (!any(id(i) == _staircaserPar.idList))
            ## already closed or never created, so report success
            success(i) = 1;
        else
            ## clear staircase information
            index = find(id(i) == _staircaserPar.idList);
            _staircaserPar.staircase(index:(end-1)) = _staircaserPar.staircase((index+1):end);
            _staircaserPar.staircase(end) = [];
            _staircaserPar.idList(index:(end-1)) = _staircaserPar.idList((index+1):end);
            _staircaserPar.idList(end) = [];
            success(i) = 1;
            if (_staircaserPar.debug >= 1)
                fprintf("%s: deleted staircase %d\n", mfilename, id);
            endif
        endif
    endfor
    argout = {success};
endfunction

###########################################################################
### Staircaser("StartTrial")

function argout = StaircaserStartTrial(argin)
    global _staircaserPar;
    helpText = sprintf("\n[success, value, track] = Staircaser(\"StartTrial\", id [, track]);\n\n");
    if (CheckForHelpRequest(helpText))
        argout = {[]};
        return;
    endif

    ## process input arguments
    nargs = numel(argin);
    CheckNumberOfInputArguments(nargs, 1, 1, helpText);
    ## process output arguments
    CheckNumberOfOutputArguments(_staircaserPar.nOutputArgs, 0, 3, helpText);

    id = argin{1};
    AssertValidId(id);
    index = find(id == _staircaserPar.idList);
    if (_staircaserPar.staircase(index).isDone)
        success = 1;
        label = 0;
        value = _staircaserPar.staircase(index).finalValue;
        if (_staircaserPar.debug >= 1)
            fprintf("%s: staircase %d is complete\n", mfilename, id);
        endif
    elseif (_staircaserPar.staircase(index).inTrial)
        success = 0;
        label = [];
        value = [];
        if (_staircaserPar.debug >= 2)
            fprintf("%s: trial already started for staircase %d\n", ...
                    mfilename, id);
        endif
    else
        _staircaserPar.staircase(index).inTrial = 1;
        ## pick the next remaining track
        currentTrack = mod(_staircaserPar.staircase(index).lastTrack, ...
                           _staircaserPar.staircase(index).nTracksRemaining) + 1;
        _staircaserPar.staircase(index).currentTrack = currentTrack;
        value = _staircaserPar.staircase(index).tracks(currentTrack).value;
        label = _staircaserPar.staircase(index).tracks(currentTrack).label;
        success = 1;
        if (_staircaserPar.debug >= 2)
            fprintf("%s: started trial for staircase %d, track %d\n", ...
                    mfilename, id, label);
        endif
    endif
    argout = {success, value, label};
endfunction # StartTrial


###########################################################################
### Staircaser("EndTrial")

function argout = StaircaserEndTrial (argin)
    global _staircaserPar;
    helpText = sprintf(["\n[success, isDone, reversal] = Staircaser(\"EndTrial\", ", ...
                        "id, response);\n\n"]);
    if (CheckForHelpRequest(helpText))
        argout = {[]};
        return;
    endif

    nargs = numel(argin);
    CheckNumberOfInputArguments(nargs, 2, 2, helpText);
    id = argin{1};
    AssertValidId(id);
    index = find(id == _staircaserPar.idList);
    response = argin{2};
    if (_staircaserPar.staircase(index).isDone)
        success = 1;
        isDone = 1;
        reversal = 0;
    elseif (!(_staircaserPar.staircase(index).inTrial))
        success = 0;
        isDone = _staircaserPar.staircase(index).isDone;
        reversal = 0;
    elseif (response < 0 || response > numel(_staircaserPar.staircase(index).steps))
        error(["invalid response code; for staircase %d, must be in range ", ...
               "[0, %d]"], id, numel(_staircaserPar.staircase(index).steps));
    else
        [isDone, reversal] = UpdateStaircase(index, response);
        success = 1;
    endif
    argout = {success, isDone, reversal};
endfunction # EndTrial

###########################################################################
### Staircaser("FinalValue")

function argout = StaircaserFinalValue (argin)
    global _staircaserPar;
    helpText = sprintf("\nvalue = Staircaser(\"FinalValue\", id);\n\n");
    if (CheckForHelpRequest(helpText));
        argout({});
        return;
    endif
    nargs = numel(argin);
    CheckNumberOfInputArguments(nargs, 1, 1, helpText);
    id = argin{1};
    AssertValidId(id);
    index = find(id == _staircaserPar.idList);
    if (isempty(_staircaserPar.staircase(index).finalValue))
        value = ComputeFinalValue(index);
    else
        value = _staircaserPar.staircase(index).finalValue;
    endif
    argout = {value};
endfunction


###########################################################################
### Staircaser("GetReversals")

function argout = StaircaserGetReversals(argin)
    global _staircaserPar;
    helpText = sprintf("\nreversals = Staircaser(\"GetReversals\", id);\n\n");
    if (CheckForHelpRequest(helpText));
        argout({});
        return;
    endif
    nargs = numel(argin);
    CheckNumberOfInputArguments(nargs, 1, 1, helpText);
    id = argin{1};
    AssertValidId(id);
    index = find(id == _staircaserPar.idList);
    if (isempty(_staircaserPar.staircase(index).reversals))
        reversals = GatherReversals(index);
    else
        reversals = _staircaserPar.staircase(index).reversals;
    endif
    argout = {reversals};
endfunction


###########################################################################
### Staircaser("Plot")

function argout = StaircaserPlot(argin)
    global _staircaserPar;
    if (_staircaserPar.debug >= 1)
        printf("Staircaser(\"Plot\") is not implemented\n");
    endif
    argout = {};
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("Progress")

function argout = StaircaserProgress(argin)

    argout = {[]};

    function Help
        fprintf(["\n[progress, stepsize] = " ...
                 "Staircaser(\"Progress\", id);\n\n"]);
    end
    if printHelp
        Help;
        return
    end

    % process input args
    nargs = numel(argin);
    if nargs < 1
        Help;
        error("not enough input arguments");
    elseif nargs > 1
        Help;
        error("too many input arguments");
    end
    id = argin{1};

    % error checking
    AssertValidId(id);

    progress = staircase(id).progress ./ staircase(id).progressTotal;
    stepsize = 1 ./ staircase(id).progressTotal;
    argout = {progress, stepsize};

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("List")

function argout = StaircaserList(argin)

function Help
fprintf(["\nid = Staircaser(\"List\");\n\n"]);
end

if printHelp
    Help;
    argout = {[]};
    return
end

% process input arguments
nargs = numel(argin);
if nargs > 0
    Help;
    error("too many input arguments");
end
% process output arguments
if nOutputArgs > 1
    Help;
    error("too many output arguments");
end

if isempty(idList) || all(idList == 0)
    argout = {[]};
else
    argout = {find(idList == 1)};
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Private support functions

function [isDone, reversal] = UpdateStaircase (index, response)
    global _staircaserPar;
    reversal = 0;

    ## extract info from staircase
    staircase = _staircaserPar.staircase(index);
    ctrack = staircase.currentTrack;
    track = staircase.tracks(ctrack);
    curval = track.value;
    lastStepDir = track.lastStepDir;
    label = track.label;
    nTrialsBase = _staircaserPar.nTrialsBase;

    ## reset staircase
    staircase.inTrial = 0;
    ## expand trial vectors, if necessary
    n = track.nTrials + 1;
    if (n > numel(track.values))
        ## need more room in values and responses array
        track.values = [track.values, nan(nTrialsBase, 1)];
        track.responses = [track.responses, nan(nTrialsBase, 1)];
        if (_staircaserPar.debug >= 1)
            fprintf("%s: expanded data matrix to length %d\n", mfilename,
                    numel(track.values));
        endif
    endif

    ## store staircase value and response for this trial
    track.values(n) = curval;
    track.responses(n) = response;
    track.nTrials = n;

    ## update staircase only if there was a valid response
    if (response != 0)
        ## compute step and step direction
        step = staircase.steps(response);
        stepDir = sign(step); % could be zero, if step is zero
        if (_staircaserPar.debug >= 2)
            fprintf("%s: took step of %0.3f on staircase %d, track %d\n",
                    mfilename, step, _staircaserPar.idList(index), label);
        endif
        ## check for reversal
        if (!isempty(lastStepDir) && stepDir != 0 && stepDir != lastStepDir)
            if (_staircaserPar.debug >= 2)
                fprintf("%s: reversal on staircase %d, track %d\n", ...
                        mfilename, _staircaserPar.idList(index), label);
            endif
            reversal = track.counter + 1;
            track.reversals(reversal) = curval;
            track.counter = reversal;
            staircase.progress = staircase.progress + 1;
        end
        ## update staircase value
        if (stepDir != 0)
            newval = curval + step;
            if (!isempty(staircase.range))
                range = staircase.range;
                if (newval > max(range))
                    newval = max(range);
                endif
                if (newval < min(range))
                    newval = min(range);
                endif
            endif
            track.lastStepDir = stepDir;
            track.value = newval;
        endif
        staircase.lastTrack = ctrack;
    endif

    ## restore track to staircase variable
    staircase.tracks(ctrack) = track;

    ## clean out completed tracks
    if track.counter >= staircase.nReversals
        ## this track is done
        if (_staircaserPar.debug >= 1)
            fprintf("%s: completed track %d for staircase %d\n", ...
                    mfilename, label, _staircaserPar.idList(index));
        endif
        nTracksRemaining = staircase.nTracksRemaining;
        if (ctrack != nTracksRemaining)
            ## move this track to the end of the list so the active ones
            ## are at the beginning
            tmp = staircase.tracks(nTracksRemaining);
            staircase.tracks(nTracksRemaining) = staircase.tracks(ctrack);
            staircase.tracks(ctrack) = tmp;
        endif
        staircase.nTracksRemaining = nTracksRemaining - 1;
        staircase.lastTrack = 0;
    endif
    if (_staircaserPar.debug >= 2)
        fprintf("%s: ended trial for staircase %d, track %d\n", ...
                mfilename, _staircaserPar.idList(index), label);
    endif
    staircase.currentTrack = 0;
    isDone = (staircase.nTracksRemaining == 0);

    ## restore staircase to global variable
    _staircaserPar.staircase(index) = staircase;

    ## close out staircase if done
    if (isDone)
        FinalizeStaircase(index);
    endif
    success = 1;
endfunction


function FinalizeStaircase (index)
    global _staircaserPar;
    _staircaserPar.staircase(index).reversals = GatherReversals(index);
    _staircaserPar.staircase(index).values = GatherValues(index);
    _staircaserPar.staircase(index).responses = GatherResponses(index);
    _staircaserPar.staircase(index).finalValue = ComputeFinalValue(index);
    _staircaserPar.staircase(index).isDone = 1;
    _staircaserPar.staircase(index).tracks = [];
endfunction


function rev = GatherReversals (index)
    global _staircaserPar;
    rev = zeros(_staircaserPar.staircase(index).nReversals,
                _staircaserPar.staircase(index).nTracks);
    for (t = 1:size(rev, 2))
        rev(:, t) = _staircaserPar.staircase(index).tracks(t).reversals;
    endfor
endfunction


function val = GatherValues (index)
    global _staircaserPar;
    nTracks = _staircaserPar.staircase(index).nTracks;
    nTrials = zeros(nTracks, 1);
    for (t = 1:nTracks)
        nTrials(t) = _staircaserPar.staircase(index).tracks(t).nTrials;
    endfor
    val = nan(max(nTrials), nTracks);
    for (t = 1:nTracks)
        if (nTrials(t) > 0)
            val(1:nTrials(t), t) = ...
                _staircaserPar.staircase(index).tracks(t).values(1:nTrials(t));
        endif
    endfor
endfunction


function resp = GatherResponses (index)
    global _staircaserPar;
    nTracks = _staircaserPar.staircase(index).nTracks;
    nTrials = zeros(nTracks, 1);
    for (t = 1:nTracks)
        nTrials(t) = _staircaserPar.staircase(index).tracks(t).nTrials;
    endfor
    resp = nan(max(nTrials), nTracks);
    for (t = 1:nTracks)
        if nTrials(t) > 0
            resp(1:nTrials(t), t) = ...
                _staircaserPar.staircase(index).tracks(t).responses(1:nTrials(t));
        endif
    endfor
endfunction


function finalValue = ComputeFinalValue (index)
    global _staircaserPar;
    if (isempty(_staircaserPar.staircase(index).reversals))
        rev = GatherReversals(index);
    else
        rev = _staircaserPar.staircase(index).reversals;
    endif
    drop = _staircaserPar.staircase(index).nReversalsDropped;
    if (drop > 0)
        rev = rev(drop+1:end, :);
    endif
    rev = reshape(rev, [numel(rev), 1]);
    if (all(isnan(rev)))
        finalValue = [];
    else
        finalValue = mean(rev(!isnan(rev)));
    endif
endfunction


function AssertValidId (id)
    global _staircaserPar;
    if (!isnumeric(id) || !any(id == _staircaserPar.idList))
        error("invalid staircase id");
    endif
endfunction


### Local Variables:
### mode:Octave
### End:
