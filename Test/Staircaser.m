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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("Delete")

function argout = StaircaserDelete(argin)

function Help
fprintf("\nsuccess = Staircaser(\"Delete\", id);\n\n");
end

if printHelp
    Help;
    argout = {[]};
end

% process input arguments
nargs = numel(argin);
if nargs < 1
    Help;
    error("not enough input arguments");
elseif nargs > 1
    Help;
    error("too many input arguments");
end
id = argin{1};
success = zeros(1, numel(id));
for i = 1:numel(id)
    if ~isnumeric(id(i)) || id(i) < 1 || id(i) > numel(idList)
        error("invalid staircase id given");
    elseif idList(id(i)) == 0;
        success(i) = 1;
    else
        % clear staircase information
        staircase(id(i)).type = [];
        staircase(id(i)).nReversals = [];
        staircase(id(i)).steps = [];
        staircase(id(i)).nReversalsDropped = [];
        staircase(id(i)).nTracks = [];
        staircase(id(i)).nTracksRemaining = [];
        staircase(id(i)).range = [];
        staircase(id(i)).inTrial = [];
        staircase(id(i)).currentTrack = [];
        staircase(id(i)).tracks = [];
        idList(id(i)) = 0;
        success(i) = 1;
    end
end
% process output arguments
if nOutputArgs > 1
    Help;
    error("too many output arguments");
end

% manage list of staircase ids
if numel(idList) >= idListBase && sum(idList) <= numel(idList) / 2.0 && ...
        all(idList(end-idListBase+1:end) == 0)
    % idList is mostly empty and has blank slots at the end, so clean it up
    i = numel(idList);
    while idList(i) == 0, i = i - 1; end
    idList = idList(1:i);
    if debug >= 1
        fprintf("%s: shrunk idList to length %d\n", mfilename, ...
                numel(idList));
    end
end

if debug >= 1
    fprintf("%s: deleted staircase %d\n", mfilename, id);
end
argout = {success};

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("StartTrial")

function argout = StaircaserStartTrial(argin)

function Help
fprintf("\n[success, value, track] = Staircaser(\"StartTrial\", id);\n\n");
end
if printHelp
    Help;
    argout = {[]};
    return
end

% process input arguments
nargs = numel(argin);
if nargs < 1
    Help;
    error("not enough input arguments");
elseif nargs > 1
    Help;
    error("too many input arguments");
end

% process output arguments
if nOutputArgs > 3
    Help;
    error("too many output arguments");
end
id = argin{1};
AssertValidId(id);
if staircase(id).isDone
    success = 1;
    label = 0;
    value = staircase(id).finalValue;
    if debug >= 1
        fprintf("%s: staircase %d is complete\n", mfilename, id);
    end
elseif staircase(id).inTrial
    success = 0;
    label = [];
    value = [];
    if debug >= 2
        fprintf("%s: trial already started for staircase %d\n", ...
                mfilename, id);
    end
else
    staircase(id).inTrial = 1;
    % pick the next remaining track
    currentTrack = mod(staircase(id).lastTrack, ...
                       staircase(id).nTracksRemaining) + 1;
    staircase(id).currentTrack = currentTrack;
    value = staircase(id).tracks(currentTrack).value;
    success = 1;
    label = staircase(id).tracks(currentTrack).label;
    if debug >= 2
        fprintf("%s: started trial for staircase %d, track %d\n", ...
                mfilename, id, label);
    end
end

argout = {success, value, label};

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("EndTrial")

function argout = StaircaserEndTrial(argin);

function Help
fprintf(["\n[success, isDone, reversal] = Staircaser(\"EndTrial\", " ...
         "id, response);\n\n"]);
end

if printHelp
    Help;
    argout = {[]};
    return
end

nargs = numel(argin);
if nargs < 2
    Help;
    error("not enough input arguments");
elseif nargs > 2
    Help;
    error("too many input arguments");
end
id = argin{1};
response = argin{2};
AssertValidId(id);
if staircase(id).isDone
    success = 1;
    isDone = 1;
    reversal = 0;
elseif ~staircase(id).inTrial
    success = 0;
    isDone = staircase(id).isDone;
    reversal = 0;
elseif response < 0 || response > numel(staircase(id).steps)
    error(["invalid response code; for staircase %d, must be in range " ...
           "[0, %d]"], id, numel(staircase(id).steps));
else
    reversal = 0;
    % reset staircase
    staircase(id).inTrial = 0;
    % extract info from staircase (to reduce repeated struct lookups)
    ctrack = staircase(id).currentTrack;
    curval = staircase(id).tracks(ctrack).value;
    lastStepDir = staircase(id).tracks(ctrack).lastStepDir;
    label = staircase(id).tracks(ctrack).label;

    % store staircase value and response for this trial
    n = staircase(id).tracks(ctrack).nTrials + 1;
    if n > numel(staircase(id).tracks(ctrack).values)
        % need more room in values and responses array
        x = nan(n + nTrialsBase - 1, 1);
        x(1:numel(staircase(id).tracks(ctrack).values)) = ...
            staircase(id).tracks(ctrack).values;
        staircase(id).tracks(ctrack).values = x;
        x = nan(n + nTrialsBase - 1, 1);
        x(1:numel(staircase(id).tracks(ctrack).responses)) = ...
            staircase(id).tracks(ctrack).responses;
        staircase(id).tracks(ctrack).responses = x;
        if debug >= 1
            fprintf("%s: expanded data matrix to length %d\n", ...
                    mfilename, numel(x));
        end
    end
    staircase(id).tracks(ctrack).values(n) = curval;
    staircase(id).tracks(ctrack).responses(n) = response;
    staircase(id).tracks(ctrack).nTrials = n;

    % update staircase only if there was a valid response
    if response ~= 0
        % compute step and step direction
        step = staircase(id).steps(response);
        stepDir = sign(step); % could be zero, if step is zero
        if debug >= 2
            fprintf("%s: took step of %0.3f on staircase %d, track %d\n", ...
                    mfilename, step, id, label);
        end
        % check for reversal
        if ~isempty(lastStepDir) && stepDir ~= 0 && stepDir ~= lastStepDir
            if debug >= 2
                fprintf("%s: reversal on staircase %d, track %d\n", ...
                        mfilename, id, label);
            end
            reversal = staircase(id).tracks(ctrack).counter + 1;
            staircase(id).tracks(ctrack).reversals(reversal) = curval;
            staircase(id).tracks(ctrack).counter = reversal;
            staircase(id).progress = staircase(id).progress + 1;
        end
        % update staircase value
        if stepDir ~= 0
            newval = curval + step;
            if ~isempty(staircase(id).range)
                range = staircase(id).range;
                if newval > max(range), newval = max(range); end
                if newval < min(range), newval = min(range); end
            end
            staircase(id).tracks(ctrack).lastStepDir = stepDir;
            staircase(id).tracks(ctrack).value = newval;
        end
        staircase(id).lastTrack = ctrack;
    end

    % clean out completed tracks
    if staircase(id).tracks(ctrack).counter >= staircase(id).nReversals
        % this track is done
        if debug >= 1
            fprintf("%s: completed track %d for staircase %d\n", ...
                    mfilename, label, id);
        end
        nTracksRemaining = staircase(id).nTracksRemaining;
        if ctrack ~= nTracksRemaining
            % move this track to the end of the list so the active ones are
            % at the beginning
            tmp = staircase(id).tracks(nTracksRemaining);
            staircase(id).tracks(nTracksRemaining) = ...
                staircase(id).tracks(ctrack);
            staircase(id).tracks(ctrack) = tmp;
        end
        staircase(id).nTracksRemaining = nTracksRemaining - 1;
        staircase(id).lastTrack = 0;
    end
    if debug >= 2
        fprintf("%s: ended trial for staircase %d, track %d\n", ...
                mfilename, id, label);
    end
    isDone = (staircase(id).nTracksRemaining == 0);
    if isDone
        FinalizeStaircase(id);
    end
    staircase(id).currentTrack = 0;
    success = 1;
end
argout = {success, isDone, reversal};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("FinalValue")

function argout = StaircaserFinalValue(argin)

function Help
fprintf("\nvalue = Staircaser(\"FinalValue\", id);\n\n");
end
if printHelp
    Help;
    argout = {[]};
    return
end

nargs = numel(argin);
if nargs < 1
    Help;
    error("not enough input arguments");
elseif nargs > 1
    Help;
    error("too many input arguments");
end
id = argin{1};
AssertValidId(id);

if isempty(staircase(id).finalValue)
    value = ComputeFinalValue(id);
else
    value = staircase(id).finalValue;
end
argout = {value};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("GetReversals")

function argout = StaircaserGetReversals(argin)

function Help
fprintf("\nreversals = Staircaser(\"GetReversals\", id);\n\n");
end
if printHelp
    Help;
    argout = {[]};
    return
end

nargs = numel(argin);
if nargs < 1
    Help;
    error("not enough input arguments");
elseif nargs > 1
    Help;
    error("too many input arguments");
end

id = argin{1};
AssertValidId(id);
if isempty(staircase(id).reversals)
    reversals = GatherReversals(id);
else
    reversals = staircase(id).reversals;
end
argout = {reversals};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Staircaser("Plot")

function argout = StaircaserPlot(argin)

argout = {[]};

function Help
fprintf("\nStaircaser(\"Plot\", id, title);\n\n");
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
elseif nargs > 2
    Help;
    error("too many input arguments");
end
id = argin{1};
if nargs < 2
    plotTitle = [];
else
    plotTitle = argin{2};
end

% error checking
AssertValidId(id);

% gather data
if isempty(staircase(id).values)
    y = GatherValues(id);
else
    y = staircase(id).values;
end
if isempty(staircase(id).finalValue)
    final = ComputeFinalValue(id);
else
    final = staircase(id).finalValue;
end
nTracks = size(y, 2);
nTrials = zeros(nTracks, 1);
for t = 1:nTracks
    nTrials(t) = max(find(~isnan(y(:, t))));
end
if ~any(nTrials > 0)
    warning("no data to plot");
    return;
end

% plot first track with data to plot, and establish axes, etc.
first = min(find(nTrials > 0));
figure;
plot(1:nTrials(first), y(1:nTrials(first), first), ...
     "-o", "LineWidth", 2, "MarkerFaceColor", "w");
hold all;
if isempty(staircase(id).range)
    axis([1, max(nTrials)+1, min(min(y)), max(max(y))]);
else
    axis([1, max(nTrials)+1, ...
          min(staircase(id).range), max(staircase(id).range)]);
end        
if ~isempty(plotTitle), title(plotTitle); end

% plot any remaining tracks that have data
offset = 0;
for t = first+1:nTracks
    if nTrials(t) == 0, continue; end
    offset = offset + .1;
    plot(offset + (1:nTrials(t)), y(1:nTrials(t), t), ...
         "-o", "LineWidth", 2, "MarkerFaceColor", "w");
end

% plot final value
plot([1, max(nTrials)+1], [final, final], "k-", "LineWidth", 3);
hold off;

end


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

function FinalizeStaircase (id)

    staircase(id).reversals = GatherReversals(id);
    staircase(id).values = GatherValues(id);
    staircase(id).responses = GatherResponses(id);
    staircase(id).finalValue = ComputeFinalValue(id);
    staircase(id).isDone = 1;
    staircase(id).tracks = [];

end


function rev = GatherReversals (id)

    AssertValidId(id);

    rev = zeros(staircase(id).nReversals, staircase(id).nTracks);
    for t = 1:size(rev, 2)
        rev(:, t) = staircase(id).tracks(t).reversals;
    end

end


function val = GatherValues (id)

    AssertValidId(id);

    nTracks = staircase(id).nTracks;
    nTrials = zeros(nTracks, 1);
    for t = 1:nTracks
        nTrials(t) = staircase(id).tracks(t).nTrials;
    end
    val = nan(max(nTrials), nTracks);
    for t = 1:nTracks
        if nTrials(t) > 0
            val(1:nTrials(t), t) = ...
                staircase(id).tracks(t).values(1:nTrials(t));
        end
    end

end


function resp = GatherResponses (id)

    AssertValidId(id);

    nTracks = staircase(id).nTracks;
    nTrials = zeros(nTracks, 1);
    for t = 1:nTracks
        nTrials(t) = staircase(id).tracks(t).nTrials;
    end
    resp = nan(max(nTrials), nTracks);
    for t = 1:nTracks
        if nTrials(t) > 0
            resp(1:nTrials(t), t) = ...
                staircase(id).tracks(t).responses(1:nTrials(t));
        end
    end

end


function finalValue = ComputeFinalValue (id)

    AssertValidId(id);

    if isempty(staircase(id).reversals)
        rev = GatherReversals(id);
    else
        rev = staircase(id).reversals;
    end

    drop = staircase(id).nReversalsDropped;
    if drop > 0
        rev = rev(drop+1:end, :);
    end
    rev = reshape(rev, [numel(rev), 1]);
    if all(isnan(rev))
        finalValue = [];
    else
        finalValue = mean(rev(~isnan(rev)));
    end
    
end


function AssertValidId (id)
    if ~isnumeric(id) || id < 1 || id > numel(idList) || idList(id) == 0
        error("invalid staircase id");
    end
end


### Local Variables:
### mode:Octave
### End:
