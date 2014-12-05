function SpeedSearch (varargin)
    global Experiment = 'test';
    global Version = '0.01';
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
    # Animation speed
    nRefreshesPerFrame = 2;
    # Size of the underlying sinewave patch
    par.gaborSize = 200;
    par.gaborRect = [0 0 par.gaborSize par.gaborSize];
    # Gabor drift speed
    phaseStep = 30;
    # Starting phase
    phase = 0;
    # Gabor frequency (between about .05 and .2 is reasonable)
    freq = .08;
    # Size of gaussian envelope
    spatialconstant = 20;
    # Sorta like contrast, but not exactly
    contrast = 100;
    # Ignored unless a parameter is set in the gabor code
    aspectratio = 1.0;
    # Angle in degrees
    angle = 0;
    gabortex = CreateProceduralGabor(par.mainWindow, par.gaborSize, par.gaborSize);
    KbReleaseWait();
    Screen('FillRect', par.mainWindow, 128);
    t = Screen('Flip', par.mainWindow);
    tNext = t + nRefreshesPerFrame * par.refreshDuration - par.slackDuration;
    phase = phase - phaseStep;
    while (1)
        phase = mod(phase + phaseStep, 360);
        Screen('FillRect', par.mainWindow, 128);
        Screen('DrawTexture', par.mainWindow, gabortex, [], [], angle, [], [], [], [], kPsychDontDoRotation, ...
               [phase+180, freq, spatialconstant, contrast]);
        t = Screen('Flip', par.mainWindow, tNext);
        tNext = t + nRefreshesPerFrame * par.refreshDuration - par.slackDuration;
        if (KbCheck)
           break;
        endif
    endwhile
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
