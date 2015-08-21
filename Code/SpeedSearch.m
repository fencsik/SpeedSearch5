function SpeedSearch (varargin)
    global Experiment = 'test';
    global Version = '0.02';
    global TestFlag = 1;
    if (nargin > 0)
        HandleInputArguments(varargin{:})
        return;
    endif
    try
        RunExperiment();
    catch
        ple();
    end_try_catch
    try
        SaveBlockData();
    catch
        ple();
    end_try_catch
    Shutdown();
    clear -global;
    clear -all;
endfunction

function RunExperiment
    Initialize();
    RunBlock();
    Deinitialize();
endfunction


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
    %% Gabor drift speed
    randomSign = randi(2, 1, nGabors);
    randomSign(randomSign == 2) = -1;
    phaseStep = randomSign .* randi([15 30], 1, nGabors);
    %% Starting phase
    phase = 0;
    %% Gabor frequency (between about .05 and .2 is reasonable)
    freq = .08;
    %% Size of gaussian envelope
    spatialconstant = 20;
    %% Sorta like contrast, but not exactly
    contrast = 100;
    %% Ignored unless a parameter is set in the gabor code
    aspectratio = 1.0;
    %% Angle in degrees
    angle = 0;
    %% Timing checks
    tLastFlip = NA;
    maxFrameDur = -1;
    sumFrameDur = 0;
    maxPrepDur = -1;
    sumPrepDur = 0;
    nFrames = 0;
    parameters = repmat([phase + 180, freq, spatialconstant, contrast]', 1, nGabors);
    KbReleaseWait();
    Screen('FillRect', par.mainWindow, 128);
    t = Screen('Flip', par.mainWindow);
    tNext = t + par.frameDuration;
    phase = phase - phaseStep;
    while (1)
        tPrepStart = GetSecs();
        parameters(1, :) = mod(parameters(1, :) + phaseStep, 360);
        Screen('FillRect', par.mainWindow, 128);
        Screen('DrawTextures', par.mainWindow, par.gabortex,
               [], par.destRect, angle, [], [], [], [], kPsychDontDoRotation, parameters);
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
            endif
            if (prepDur > maxPrepDur)
                maxPrepDur = prepDur;
            endif
        endif
        tLastFlip = t;
        if (KbCheck)
           break;
        endif
    endwhile
    if (nFrames > 0)
        fprintf('%0.0f frames completed\n', nFrames);
        fprintf('Average frame duration = %0.6f ms\n', 1000 * sumFrameDur / nFrames);
        fprintf('Maximum frame duration = %0.6f ms\n', 1000 * maxFrameDur);
        fprintf('Average prep duration  = %0.6f ms\n', 1000 * sumPrepDur / nFrames);
        fprintf('Maximum prep duration  = %0.6f ms\n', 1000 * maxPrepDur);
    endif
endfunction

function SaveBlockData ()
endfunction



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization and Shutdown Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Initialize ()
    InitializePreGraphics();
    InitializeGraphics();
    InitializePostGraphics();
endfunction

function InitializePreGraphics ()
    AssertOpenGL();
    KbName('UnifyKeyNames');
    global par = struct();

    %% Stimulus layout
    par.nClusters = 3;
    par.nStimuliPerCluster = 4;
    par.clusterSpacingDenominator = 16; % >= par.nClusters * par.nStimuliPerCluster
    par.displayRadius = 360;

    %% Animation speed
    par.nRefreshesPerFrame = 1;

    %% Size of the gabor patch
    par.gaborSize = 200;
    par.gaborRect = [0 0 par.gaborSize par.gaborSize];
endfunction

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
endfunction

function InitializePostGraphics ()
    global par
    Screen('BlendFunction', par.mainWindow, GL_ONE, GL_ONE);
    %%Screen('BlendFunction', par.mainWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    HideCursor();

    [par.screenCenterX, par.screenCenterY] = RectCenter(par.mainWindowRect);

    %% Create gabor texture
    par.gabortex = CreateProceduralGabor(par.mainWindow, par.gaborSize, par.gaborSize);

    %% calculate frame durations and number of frames
    par.refreshDuration = Screen('GetFlipInterval', par.mainWindow);
    par.slackDuration = par.refreshDuration / 2.0;
    par.frameDuration = par.nRefreshesPerFrame * par.refreshDuration - par.slackDuration;
endfunction

function Deinitialize ()
end

function Shutdown ()
    Priority(0);
    fclose('all');
    ShutdownGraphics();
endfunction

function ShutdownGraphics ()
    ShowCursor();
    Screen('CloseAll');
endfunction


function HandleInputArguments (varargin)
    global InputArguments
    if (nargin > 0)
        InputArguments.nGabors = varargin{1};
    endif
endfunction


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
