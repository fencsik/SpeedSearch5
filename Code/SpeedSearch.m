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


###########################################################################
### Block-Level Functions
###########################################################################

function RunBlock ()
    global par
    ## Animation speed
    nRefreshesPerFrame = 2;
    ## Size of the underlying sinewave patch
    par.gaborSize = 200;
    par.gaborRect = [0 0 par.gaborSize par.gaborSize];
    ## Define two destination rects
    centeredGaborRect = CenterRect(par.gaborRect, par.mainWindowRect);
    par.destRect = [OffsetRect(centeredGaborRect, -1 * par.gaborSize, 0)', ...
                    OffsetRect(centeredGaborRect, 1 * par.gaborSize, 0)', ...
                    OffsetRect(centeredGaborRect, 0, -1 * par.gaborSize)', ...
                    OffsetRect(centeredGaborRect, 0, 1 * par.gaborSize)'];
    ## Gabor drift speed
    phaseStep = [25, -25, 50, -50];
    ## Starting phase
    phase = 0;
    ## Gabor frequency (between about .05 and .2 is reasonable)
    freq = .08;
    ## Size of gaussian envelope
    spatialconstant = 20;
    ## Sorta like contrast, but not exactly
    contrast = 100;
    ## Ignored unless a parameter is set in the gabor code
    aspectratio = 1.0;
    ## Angle in degrees
    angle = 0;
    ## Timing checks
    tLastFlip = NA;
    maxFrameDur = -1;
    sumFrameDur = 0;
    maxPrepDur = -1;
    nFrames = 0;
    parameters = repmat([phase + 180, freq, spatialconstant, contrast]', 1, 4);
    gabortex = CreateProceduralGabor(par.mainWindow, par.gaborSize, par.gaborSize);
    KbReleaseWait();
    Screen('FillRect', par.mainWindow, 128);
    t = Screen('Flip', par.mainWindow);
    tNext = t + nRefreshesPerFrame * par.refreshDuration - par.slackDuration;
    phase = phase - phaseStep;
    while (1)
        tPrepStart = GetSecs();
        parameters(1, :) = mod(parameters(1, :) + phaseStep, 360);
        Screen('FillRect', par.mainWindow, 128);
        Screen('DrawTextures', par.mainWindow, [gabortex, gabortex, gabortex, gabortex],
               [], par.destRect, angle, [], [], [], [], kPsychDontDoRotation, parameters);
        tPrepEnd = GetSecs();
        t = Screen('Flip', par.mainWindow, tNext);
        tNext = t + nRefreshesPerFrame * par.refreshDuration - par.slackDuration;
        if (!isna(tLastFlip))
            dur = t - tLastFlip;
            sumFrameDur = sumFrameDur + dur;
            nFrames = nFrames + 1;
            if (dur > maxFrameDur)
                maxFrameDur = dur;
            endif
            dur = tPrepEnd - tPrepStart;
            if (dur > maxPrepDur)
                maxPrepDur = dur;
            endif
        endif
        tLastFlip = t;
        if (KbCheck)
           break;
        endif
    endwhile
    if (nFrames > 0)
        printf("%0.0f frames completed\n", nFrames);
        printf("Average frame duration = %0.6f ms\n", 1000 * sumFrameDur / nFrames);
        printf("Maximum frame duration = %0.6f ms\n", 1000 * maxFrameDur);
        printf("Maximum prep duration  = %0.6f ms\n", 1000 * maxPrepDur);
    endif
endfunction

function SaveBlockData ()
endfunction



###########################################################################
### Initialization and Shutdown Functions
###########################################################################

function Initialize ()
    InitializePreGraphics();
    InitializeGraphics();
    InitializePostGraphics();
endfunction

function InitializePreGraphics ()
    AssertOpenGL();
    KbName('UnifyKeyNames');
    global par = struct();
endfunction

function InitializeGraphics ()
    global par;
    screenNumber=max(Screen('Screens'));

    ## # Old style of opening windows
    ## [par.mainWindow, par.mainWindowRect] = ...
    ##     Screen('OpenWindow', screenNumber, par.backgroundColor, [], 32, 2);
    ## Screen('BlendFunction', par.mainWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    # New style of opening windows
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    [par.mainWindow, par.mainWindowRect] = PsychImaging('OpenWindow', screenNumber, 128);
endfunction

function InitializePostGraphics ()
    global par
    Screen('BlendFunction', par.mainWindow, GL_ONE, GL_ONE);
    #Screen('BlendFunction', par.mainWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    HideCursor();

    % calculate frame durations and number of frames
    par.refreshDuration = Screen('GetFlipInterval', par.mainWindow);
    par.slackDuration = par.refreshDuration / 2.0;
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


### Local Variables:
### mode:Octave
### End:
