function [] = BackCorr_after_Analyze_ALL_channels(posI, TFperNom, myChannel,zSlices,zfocus,psfI)


    %% Declare Variables
    % setup path for blank images
    blankPath = [];
        
    loc_path = pwd;
    all_nm = [647 561 488];
    numRnds = 81;
    channel = myChannel + 1;
    channel_nm = all_nm(channel);
    my_pos = posI-1;
    
    
    %% Grab FISH images
    [Image_647] = ...
        grab_ims(loc_path, my_pos, numRnds, myChannel, zSlices);
    

    %% Grab DAPI images
    DAPI_Channel = 3;
    [Image_DAPI] = ...
        grab_ims(loc_path, my_pos, numRnds, DAPI_Channel, zSlices);
    
    
    %% use DAPI to find approximate tform between images...
    fprintf('DAPI aligning all hybs and Z-slices... \r')
    bead_im =Image_DAPI{1}(:,:,zfocus);
    trans = zeros(numRnds,2);
    tic
    for i = 1:numRnds
        my_im = Image_DAPI{i}(:,:,zfocus);
        TFORM_final(i) = imregcorr(my_im,bead_im,'translation');
        trans(i,:) = TFORM_final(i).T(3,1:2);
    end
    toc
    clear Image_DAPI
    
    
    %% check bead bleaching, Z-focus values, and tform estimates
    figure
    plot(trans(:,1),'ro')
    hold on
    plot(trans(:,2),'mo')
    legend('x est','y est','Location','best')
    savefig([loc_path sprintf('\\Analysis_Details_NO_FISH_RCE\\Bead_Alignment\\Pos%.0f_%.0fnm.fig',my_pos, channel_nm)])
    close all
    
    
    %% clean up images
    fprintf('Cleaning up images through Miji...\n')
    parfor hyb = 1:81
        tic
        hyb
        locIm = Image_647{hyb};
        myIm = locIm(:,:,zfocus);
        [Image_647_clean{hyb}] = deconvlucy(myIm,psfI);
    end
    save([ loc_path '\Analysis_Details_NO_FISH_RCE\Back_Subtract_Z_focused\' sprintf('Deconvolved_Pos%.0f_%.0fnm.mat',my_pos, channel_nm)],'Image_647_clean')

    
    %% Correct for laser
    % load blank image
    fprintf('Correcting for laser...\n')
    myChannel = channel-1;
    [Blank_Image_647] = ...
        grab_ims_Blanks(blankPath, 0:6, myChannel, zSlices);

    corrections = back_correct_multi_pos(Blank_Image_647);
    tic
    parfor hyb = 1:81
        hyb
            locIm = Image_647_clean{hyb};
            focusBeadImage{hyb} = locIm;
            BC_im{hyb} = double(locIm)./corrections;
    end
    toc
    tic
    for hyb = 1:80:81
        fprintf(sprintf('Save backcorrected %.0f ...\n',hyb)) 
        myIm = BC_im{hyb};
        save([loc_path filesep 'Analysis_Details' filesep 'Back_Subtract_Z_focused' filesep...
            sprintf('Back_Corrected_Pos%.0f_hyb%.0f_%.0fnm.mat',my_pos, hyb, channel_nm)], 'myIm')
    end
    toc
    Image_647_clean = BC_im;
    
    
    %% find FISH dots in all rounds, using rolling ball from miji
    fprintf('Finding FISH dots in all hyb rounds...\n')
    load(sprintf('threshold_factor_%.0f.mat',channel_nm))
    
    threshold_factor = threshold_factor * TFperNom;
    HCRorFISH = 1; %FISH
    debug = 0;
    FISHdots = cell(81,1);
    Pixels = 3;
    tic

    for hyb = 1:81
        myIm = Image_647_clean{hyb};
        threshold = threshold_factor(hyb);
        [dotInt{hyb},FISHdots{hyb}] = findDotsBarcodeV2({myIm}, threshold, HCRorFISH,debug);
        srDots{hyb} = FISHdots{hyb}.channels;
    end
    toc
    
    
    %% keep only dots in segmented cells
    fprintf('Removing dots outside segmented cells ...\n')
    PathName = [loc_path filesep '..' filesep];
    fullpath = [PathName filesep 'RoiSet_Pos' num2str(my_pos) ];
    vertex = selfseg(fullpath);
    for hyb = 1:81
        loc_dots = FISHdots{hyb}.channels(:,1:2);
        keep_I = zeros(length(loc_dots),1);
        for i = 1:length(vertex)
            keep_I = keep_I + ...
                inpolygon(loc_dots(:,1),loc_dots(:,2),vertex(i).x,vertex(i).y);
        end
        keep_I = find(keep_I);
        FISHdots{hyb}.channels = FISHdots{hyb}.channels(keep_I,:);
        srDots{hyb} = FISHdots{hyb}.channels;
    end
    
    
    %% super resolve with radial center applied to raw image all dots
    fprintf('SR all dots with RCE ...\n')
    Pixels = 3;
    tic
    for hyb = 1:81
        locDots = FISHdots{hyb}.channels;
        locIm = Image_647_clean{hyb};
        xc = [];
        yc = [];
        for dotID = 1:length(locDots)
            I = double((locIm((-Pixels:Pixels)+locDots(dotID,2),(-Pixels:Pixels)+locDots(dotID,1))));
            [xc(dotID),yc(dotID)] = radialcenter(I);
        end
        xyc = locDots(:,1:2) + [xc; yc]' - [Pixels, Pixels]-1;
        numdots(hyb) = length(xyc);
        srDots{hyb} = xyc;
    end
    
  
    %% build hyb labelled dot lists
    tic
    fprintf('Build hyb labelled dot lists...\n')
    dotsWid = {}; % dots with ID of actual hyb (i.e. image index -1)
    for bcI = 1:4
        myDots = [];
        for hybI = 1:20
            gHybI = (bcI-1)*20+hybI;
            hDots = srDots{gHybI};
            locDots = transformPointsForward( TFORM_final(gHybI),hDots(:,1:2));
            myDots = [myDots; locDots ones(size(locDots,1),1)*(hybI)];
        end
        dotsWid(bcI) = {myDots};
    end
    FISH_only = dotsWid;
    
    
    %% Save FISH dots
	FISH_only
    fprintf('saving results...\n')
    save([loc_path filesep sprintf('%.0f_Analysis',channel_nm) filesep sprintf('Analysis_Details_NO_FISH_RCE_%.1f',TFperNom) filesep 'extractedData' filesep...
        sprintf('FISH_only_Pos%.0f_%.0fnm_results.mat',my_pos,channel_nm)],'FISH_only')

    fprintf('results saved...\n')
% end
