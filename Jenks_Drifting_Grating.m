function angle_list = Jenks_Drifting_Grating(start_angle, cyclespersecond, f, drawmask, gratingsize,repetition,duration)
%hacked together program to do drifting grating. 
%Jenks_Drifting_Grating(start_angle, cyclespersecond, f, drawmask, gratingsize,repetition,duration)
%start_angle is a vector of angles you wish to present
%cyclespersecond is the speed of the drifting grating
%f is the spatial frequency in cycles per pixel
%drawmask can fit a gaussian mask (set to 1 if you want that)
%grating size is the size of the total stim screen
%repetition is the number of times you want to repeat each stimuli
%duration is stimulus duration 
%randomized stimuli order for each rep is shown in command line
%can set variable output to a matrix of the stimuli order
%(repetitionxangle)
%
%Will send stim timing information tthrough USB device
%USB_1208FS_Plus_Interface if connected.
%
% Made by Kyle Jenks for the use of the Shepherd lab at the University of Utah, 10-28-15


% Setup IO hardware interface to send data to Plexon
plxInterface = plxHWInterface;

% Tell Plexon to start recording
    plxInterface.startRecording;
    
try
%% Set defaults
if nargin < 5
    gratingsize = [];
end
if isempty(gratingsize)
    % By default the visible grating is 400 pixels by 400 pixels in size:
    gratingsize = 4000;
end
if nargin < 4
    drawmask = [];
end
if isempty(drawmask)
    % By default, we do not mask the grating by a gaussian transparency mask: 1
    % fits a gaussian mask
    drawmask=0;
end;
if nargin < 3
    f = [];
end
if isempty(f)
    % Grating cycles/pixel: By default 0.005 cycles per pixel.
    f=0.005;
end;
if nargin < 2
    cyclespersecond = [];
end
if isempty(cyclespersecond)
    % Speed of grating in cycles per second: 1 cycle per second by default.
    cyclespersecond=2;
end;
if nargin < 1
    start_angle = [];
end
if isempty(start_angle)
    % Angle of the grating: vector of degrees.
    start_angle=[0,30,60,90,120,150,180,210,240,270,300,330];
end;
if nargin < 1
    repetition = [];
end
if isempty(repetition)
    % repetition of each stimulus. Default 5.
    repetition=5;
end;
if nargin < 1
    duration = [];
end
if isempty(duration)
    % duration of each stimulus. Default 5 seconds.
    duration=5;
end;


%% Start psychotoolbox
% Define Half-Size of the grating image.
texsize=gratingsize / 2;
%open
AssertOpenGL;
%display on max screen
screens=Screen('Screens');
screenNumber=max(screens);
%set white and black
white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);
%set grey
gray=round((white+black)/2);
if gray == white
		gray=white / 2;
end
   
    % Contrast 'inc'rement range for given white and gray values:
	inc=white-gray;

    % Open a double buffered fullscreen window and set default background
	% color to gray:
	[w screenRect]=Screen('OpenWindow',screenNumber, gray);
    if drawmask
        % Enable alpha blending for proper combination of the gaussian aperture
        % with the drifting sine grating:
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    end
    
% Query duration of one monitor refresh interval:
    ifi=Screen('GetFlipInterval', w);
    
%% Calculate grating
% First we compute pixels per cycle, rounded up to full pixels, as we
    % need this to create a grating of proper size below:
    p=ceil(1/f);
    
    % Also need frequency in radians:
    fr=f*2*pi;
    
    % This is the visible size of the grating. It is twice the half-width
    % of the texture plus one pixel to make sure it has an odd number of
    % pixels and is therefore symmetric around the center of the texture:
    visiblesize=2*texsize+1;

    % Create one single static grating image:
    x = meshgrid(-texsize:texsize + p, 1);
    
    % Compute actual cosine grating:
    grating=gray + inc*cos(fr*x);

    % Store 1-D single row grating in texture:
    gratingtex=Screen('MakeTexture', w, grating);

    % Create a single gaussian transparency mask and store it to a texture
    mask=ones(2*texsize+1, 2*texsize+1, 2) * gray;
    [x,y]=meshgrid(-1*texsize:1*texsize,-1*texsize:1*texsize);
    mask(:, :, 2)=white * (1 - exp(-((x/90).^2)-((y/90).^2)));
    masktex=Screen('MakeTexture', w, mask);

    % Query maximum useable priorityLevel on this system:
	priorityLevel=MaxPriority(w); %#ok<NASGU>

    % We don't use Priority() in order to not accidentally overload older
    % machines that can't handle a redraw every 40 ms. If your machine is
    % fast enough, uncomment this to get more accurate timing.
    Priority(priorityLevel);
    
    % Definition of the drawn rectangle on the screen:
    % Compute it to  be the visible size of the grating, centered on the
    % screen:
    dstRect=[0 0 visiblesize visiblesize];
    dstRect=CenterRect(dstRect, screenRect);
    
    waitframes = 1;
    
    % Translate frames into seconds for screen update interval:
    waitduration = waitframes * ifi;
    
    % Recompute p, this time without the ceil() operation from above.
    % Otherwise we will get wrong drift speed due to rounding errors!
    p=1/f;  % pixels/cycle    

    % Translate requested speed of the grating (in cycles per second) into
    % a shift value in "pixels per frame", for given waitduration: This is
    % the amount of pixels to shift our srcRect "aperture" in horizontal
    % directionat each redraw:
    shiftperframe= cyclespersecond * p * waitduration;

% Generate matrix to store stim order info 
angle_list=zeros(repetition,length(start_angle));
    

    %% loop to do multiple presentations of same set of stimuli 
    for jj=1:repetition
    
    %randomize angle order
randomarray=randperm(length(start_angle));
angle=start_angle(randomarray)

%save order in output matrix
angle_list(jj,:)=angle;   

%% Start presentations 
    for ii=1:length(angle)
        pause(duration);
    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp as timing baseline for our redraw loop:
    vbl=Screen('Flip', w);
    
    %Send pulse to USB interface
    plxInterface.setEventWord(1);
    
    %define length of stimulus in seconds
    vblendtime = vbl + duration;
    i=0;
    
    % Animationloop:
    while(vbl < vblendtime)
        % Shift the grating by "shiftperframe" pixels per frame:
        % the mod'ulo operation makes sure that our "aperture" will snap
        % back to the beginning of the grating, once the border is reached.
        xoffset = mod(i*shiftperframe,p);
        i=i+1;
        % Define shifted srcRect that cuts out the properly shifted rectangular
        % area from the texture: 
        srcRect=[xoffset 0 xoffset + visiblesize visiblesize];
        % Draw grating texture, rotated by "angle":
        Screen('DrawTexture', w, gratingtex, srcRect, dstRect, angle(ii));
        if drawmask==1
            % Draw gaussian mask over grating:
            Screen('DrawTexture', w, masktex, [0 0 visiblesize visiblesize], dstRect, angle(ii));
        end;
        % Flip 'waitframes' monitor refresh intervals after last redraw.
        % Providing this 'when' timestamp allows for optimal timing
        % precision in stimulus onset, a stable animation framerate and at
        % the same time allows the built-in "skipped frames" detector to
        % work optimally and report skipped frames due to hardware
        % overload:
        vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
        % Abort  if any key is pressed:
        if KbCheck
            break;
        end;
     
    end
    
    % Now fill the screen gray, allows insterstimulus interval. 
Screen('FillRect', w, gray);
% Flip to the screen
Screen('Flip', w);
%end square wave
plxInterface.setEventWord(0);
    end
    end;
% Tell Plexon to pause recording
%set to flop
        plxInterface.setEventWord(0);
    plxInterface.stopRecording;
    
    % Restore normal priority scheduling in case something else was set
    % before:
    Priority(0);
	
	%The same commands wich close onscreen and offscreen windows also close
	%textures.
	Screen('CloseAll');

catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end %try..catch..

% We're done!
return;
