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
    [x, y] = GetUniformStimulusLocations(nGabors);
    par.destRect = CenterRectOnPoint(centeredGaborRect, x, y)';
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

    %% Stimulus layout
    par.nClusters = 3;
    par.nStimuliPerCluster = 4;
    par.clusterSpacingDenominator = 16; % >= par.nClusters * par.nStimuliPerCluster
    par.displayRadius = 300;

    %% Animation speed
    par.nRefreshesPerFrame = 1;

    %% Size of the gabor patch
    par.gaborSize = 65;
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
    freq = .05;
    %% Size of gaussian envelope
    spatialconstant = par.gaborSize / 5;
    %% Sorta like contrast, but not exactly
    contrast = 10.0;
    %% Ignored unless a parameter is set in the gabor code
    aspectratio = 1.0;
    par.gaborBaseVector = [phase; freq; spatialconstant; contrast];
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
end


%%% Local Variables:
%%% mode:Octave
%%% End:
