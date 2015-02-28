function ProceduralGaborSizing (size=1000)
    # ProceduralGaborSizing is a function to play around with the
    # relationship between a procedural gabor's parameters and it's
    # physical dimensions on the screen.
    #
    # It draws a gabor patch, with size determined by the one argument. The
    # lines drawn over the gabor are spaced 50 pixels apart for reference.
    #
    # From playing around with this function, it appears that the gabor's
    # frequency parameter is, as it claims, the number of cycles per pixel.
    # So, a frequency of .01 means a wavelength of 100 pixels.  See the
    # following site for Mario's original forum post about how this
    # function works:
    #
    # https://groups.yahoo.com/neo/groups/psychtoolbox/conversations/topics/9174
    #
    # (Note that this link works, whereas the one in the
    # CreateProceduralGabor help text doesn't.)
    try
        AssertOpenGL;
        oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
        oldSupressAllWarnings = Screen('Preference', ...
                                       'SuppressAllWarnings', 1);
        RunTest(size);
    catch
        ple();
    end_try_catch
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
endfunction

function RunTest (gaborSize)
    length = 100;
    backgroundColor = 128;
    gaborFreq = .02;
    gaborPhase = 0;
    bigStep = .010;
    smallStep = .001;
    phaseStep = 10;
    KbName('UnifyKeyNames');
    keyBigStepUp = KbName('UpArrow');
    keyBigStepDown = KbName('DownArrow');
    keySmallStepUp = KbName('RightArrow');
    keySmallStepDown = KbName('LeftArrow');
    keyPhaseUp = KbName('p');
    keyPhaseDown = KbName('o');
    keyExit = KbName('ESCAPE');
    screenNumber = max(Screen('Screens'));
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    [mainWindow, mainWindowRect] = PsychImaging('OpenWindow', screenNumber, backgroundColor);
    Screen('BlendFunction', mainWindow, GL_ONE, GL_ONE);
    gaborTex = CreateProceduralGabor(mainWindow, gaborSize, gaborSize);
    gaborAngle = 0;
    spatialConstant = gaborSize / 10;
    contrast = 200;

    [cx, cy] = RectCenter(mainWindowRect);
    mrt = mainWindowRect(RectTop);
    mrb = mainWindowRect(RectBottom);
    xy = [[cx - length; mrt], [cx - length; mrb], ...
          [cx - length / 2; mrt], [cx - length /2; mrb], ...
          [cx; mrt], [cx; mrb], ...
          [cx + length / 2; mrt], [cx + length / 2; mrb], ...
          [cx + length; mrt], [cx + length; mrb]];
    done = 0;
    while (!done)
        Screen('FillRect', mainWindow, backgroundColor);
        DrawFormattedText(mainWindow, ...
                          sprintf('frequency = %0.4f\nphase     = %0.1f degrees', ...
                                  gaborFreq, gaborPhase));
        Screen('DrawTexture', mainWindow, gaborTex, [], [], gaborAngle, ...
               [], [], [], [], kPsychDontDoRotation, ...
               [gaborPhase, gaborFreq, spatialConstant, contrast]);
        Screen('DrawLines', mainWindow, xy, 1, [0 150 0]);
        Screen('Flip', mainWindow);
        [t, keyCode] = KbPressWait();
        switch find(keyCode)(1)
            case keyExit
                done = 1;
            case keyBigStepUp
                gaborFreq = gaborFreq + bigStep;
            case keyBigStepDown
                gaborFreq = gaborFreq - bigStep;
            case keySmallStepUp
                gaborFreq = gaborFreq + smallStep;
            case keySmallStepDown
                gaborFreq = gaborFreq - smallStep;
            case keyPhaseUp
                 gaborPhase = mod(gaborPhase + phaseStep, 360);
            case keyPhaseDown
                 gaborPhase = mod(gaborPhase - phaseStep, 360);
        endswitch
    endwhile
endfunction

### Local Variables:
### mode:Octave
### End:
