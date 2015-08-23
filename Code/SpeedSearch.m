function SpeedSearch (varargin)
    global Experiment = 'test';
    global Version = '0.02';
    global TestFlag = 1;
    if (nargin > 0)
        HandleInputArguments(varargin{:})
        return;
    end
    try
        RunExperiment();
    catch
        ple();
    end
    try
        SaveBlockData();
    catch
        ple();
    end
    Shutdown();
    clear -global;
    clear -all;
end

function RunExperiment
    Initialize();
    RunBlock();
    Deinitialize();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Block-Level Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RunBlock ()
    global par
    nGabors = 12;
    %% Define destination rects
    centeredGaborRect = CenterRect(par.gaborRect, par.mainWindowRect);
    [x, y] = GetClusteredStimulusLocations(nGabors);
    par.destRect = CenterRectOnPoint(centeredGaborRect, x, y)';
    par.destRect = par.destRect(:, randperm(nGabors));
    %% Initialize gabors
    InitializeGaborsForTrial(nGabors);
    %% Gabor drift speed
    randomSign = randi(2, 1, nGabors);
    randomSign(randomSign == 2) = -1;
    phaseStep = randomSign .* randi([15 30], 1, nGabors);
    %% Timing checks
    tLastFlip = NA;
    maxFrameDur = -1;
    sumFrameDur = 0;
    maxPrepDur = -1;
    sumPrepDur = 0;
    nFrames = 0;
    KbReleaseWait();
    Screen('FillRect', par.mainWindow, 128);
    t = Screen('Flip', par.mainWindow);
    tNext = t + par.frameDuration;
    while (1)
        tPrepStart = GetSecs();
        DriftGabors(phaseStep);
        Screen('FillRect', par.mainWindow, 128);
        DrawGabors();
        tPrepEnd = GetSecs();
        t = Screen('Flip', par.mainWindow, tNext);
        tNext = t + par.frameDuration;
        if (!isna(tLastFlip))
            frameDur = t - tLastFlip;
            prepDur = tPrepEnd - tPrepStart;
            sumFrameDur = sumFrameDur + frameDur;
            sumPrepDur = sumPrepDur + prepDur;
            nFrames = nFrames + 1;
            if (frameDur > maxFrameDur)
                maxFrameDur = frameDur;
            end
            if (prepDur > maxPrepDur)
                maxPrepDur = prepDur;
            end
        end
        tLastFlip = t;
        if (KbCheck)
           break;
        end
    end
    if (nFrames > 0)
        fprintf('%0.0f frames completed\n', nFrames);
        fprintf('Average frame duration = %0.6f ms\n', 1000 * sumFrameDur / nFrames);
        fprintf('Maximum frame duration = %0.6f ms\n', 1000 * maxFrameDur);
        fprintf('Average prep duration  = %0.6f ms\n', 1000 * sumPrepDur / nFrames);
        fprintf('Maximum prep duration  = %0.6f ms\n', 1000 * maxPrepDur);
    end
end

function SaveBlockData ()
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization and Shutdown Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Initialize ()
    InitializePreGraphics();
    InitializeGraphics();
    InitializePostGraphics();
end

function InitializePreGraphics ()
    AssertOpenGL();
    KbName('UnifyKeyNames');
    global par = struct();

    %% Load settings
    LoadSettingsFile('Settings.txt');

    %% Size of the gabor patch
    par.gaborRect = [0 0 par.gaborSize par.gaborSize];
end

function InitializeGraphics ()
    global par;
    screenNumber=max(Screen('Screens'));

    %% % Old style of opening windows
    %% [par.mainWindow, par.mainWindowRect] = ...
    %%     Screen('OpenWindow', screenNumber, par.backgroundColor, [], 32, 2);
    %% Screen('BlendFunction', par.mainWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    %% New style of opening windows
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    [par.mainWindow, par.mainWindowRect] = PsychImaging('OpenWindow', screenNumber, 128);
end

function InitializePostGraphics ()
    global par
    Screen('BlendFunction', par.mainWindow, GL_ONE, GL_ONE);
    %%Screen('BlendFunction', par.mainWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    HideCursor();

    [par.screenCenterX, par.screenCenterY] = RectCenter(par.mainWindowRect);

    %% Create gabor texture
    InitializeGabors();

    %% calculate frame durations and number of frames
    par.refreshDuration = Screen('GetFlipInterval', par.mainWindow);
    par.slackDuration = par.refreshDuration / 2.0;
    par.frameDuration = par.nRefreshesPerFrame * par.refreshDuration - par.slackDuration;
end

function Deinitialize ()
end

function Shutdown ()
    Priority(0);
    fclose('all');
    ShutdownGraphics();
end

function ShutdownGraphics ()
    ShowCursor();
    Screen('CloseAll');
end


function HandleInputArguments (varargin)
    global InputArguments
    if (nargin > 0)
        InputArguments.nGabors = varargin{1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gabor Management
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitializeGabors ()
    global par
    par.gabortex = CreateProceduralGabor(par.mainWindow, par.gaborSize, par.gaborSize);

    %% Starting phase
    phase = 0;
    %% Gabor frequency (between about .05 and .2 is reasonable)
    freq = par.gaborFrequency;
    %% Size of gaussian envelope
    spatialconstant = par.gaborSize / 5;
    %% Sorta like contrast, but not exactly
    amplitudeMultiplier = par.gaborAmplitudeMultiplier;
    %% Ignored unless a parameter is set in the gabor code
    aspectratio = 1.0;
    par.gaborBaseVector = [phase; freq; spatialconstant; amplitudeMultiplier];
    par.gaborVectorPhaseIndex = 1;
end

function InitializeGaborsForTrial(nGabors)
    global par
    par.gaborTrialVector = repmat(par.gaborBaseVector, 1, nGabors);
    par.gaborTrialVector(par.gaborVectorPhaseIndex, :) = randi(360, 1, nGabors) - 1;
end

function DriftGabors(increment)
    global par
    x = par.gaborTrialVector(par.gaborVectorPhaseIndex, :);
    x = mod(x + increment, 360);
    par.gaborTrialVector(par.gaborVectorPhaseIndex, :) = x;
end

function DrawGabors()
    global par
    %% Angle in degrees
    angle = 0;
    Screen('DrawTextures', par.mainWindow, par.gabortex,
           [], par.destRect, angle, [], [], [], [], kPsychDontDoRotation,
           par.gaborTrialVector);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stimulus Drawing Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y] = GetUniformStimulusLocations (nStimuli)
    global par;
    start = 2 * pi * rand();
    spacing = 2  * pi / nStimuli;
    angles = start + (0:(nStimuli-1)) * spacing;
    x = par.screenCenterX + par.displayRadius * sin(angles');
    y = par.screenCenterY - par.displayRadius * cos(angles');
end

function [x, y] = GetClusteredStimulusLocations (nStimuli)
    global par;
    nStimulusPositions = par.nClusters * par.nStimuliPerCluster;
    if (nStimuli > nStimulusPositions)
        error('clusters cannot support more than %d stimuli (%d requested)', ...
              nStimulusPositions, nStimuli);
    end
    x = nan(nStimuli, 1);
    y = x;
    nClustersNeeded = ceil(nStimuli / par.nStimuliPerCluster);
    clusterCenter = 2 * pi * rand;
    interClusterSpacing = 2 * pi / par.nClusters;
    interStimulusSpacing = 2 * pi / par.clusterSpacingDenominator;
    counter = 1;
    for i = 1:nClustersNeeded
        theta = clusterCenter + interStimulusSpacing * ...
                ((1:par.nStimuliPerCluster) - (par.nStimuliPerCluster + 1) / 2);
        for j = 1:par.nStimuliPerCluster
            x(counter) = par.screenCenterX + par.displayRadius * sin(theta(j));
            y(counter) = par.screenCenterY - par.displayRadius * cos(theta(j));
            counter = counter + 1;
        end
        clusterCenter = clusterCenter + interClusterSpacing;
    end
    x = x(1:nStimuli);
    y = y(1:nStimuli);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Settings File Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LoadSettingsFile (settingsFileName)
    fid = fopen(settingsFileName, 'r');
    if fid != -1
        settingsArray = ProcessSettingsFile(fid);
        MakeSettingsGlobal(settingsArray);
        CloseFile(fid);
    else
        fprintf('Could not open %s\n', settingsFileName);
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

function success = CloseFile(fid)
    x = fclose(fid);
    if (x == 0)
        success = 1;
    else
        success = 0;
    end
end

function MakeSettingsGlobal (settingsArray)
    global par
    fn = fieldnames(settingsArray);
    for (i = 1:numel(fn))
        par.(fn{i}) = settingsArray.(fn{i});
    end
end


%%% Local Variables:
%%% mode:Octave
%%% End:
