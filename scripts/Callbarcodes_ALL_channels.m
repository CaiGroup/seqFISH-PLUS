function [] = Callbarcodes_ALL_channels(posI, TFperNom, radius, myChannel)

    % load barcode codebook
    barcodePath = [];
    
    %% Delcare Variables
    load(barcodePath)
    loc_pat = pwd;
    all_nm = [647 561 488];
    channel = myChannel + 1;
    channel_nm = all_nm(channel);
    FISH_TF = 1;
    bead_TF = 100;
    Path = pwd;
    numRnds = 81;
    zSlices = 5;

    % %%%%%%%%%%%%%%%%% adapt output for Sheel's code %%%%%%%%%%%%%%%%%%%%%%%%%%%
    my_pos = posI-1
    tic
    % load points
    load([loc_pat filesep sprintf('%.0f_Analysis',channel_nm)  filesep sprintf('Analysis_Details_NO_FISH_RCE_%.1f',TFperNom) filesep 'extractedData' filesep sprintf('FISH_only_Pos%.0f_%.0fnm_results.mat',my_pos, channel_nm)])

    toc
    
    hyb = 4; %bcR
    channum = 20; %hybs
    alloweddiff = 1;
    conthresh = 1;
    
    close all
    fprintf('Converting dot organization...\n')
    my_pos
    tic
    for bcR = 1:4
        locDots = FISH_only{bcR};
        for channel = 1:20
            IDs = find(locDots(:,3) == channel);
            points(bcR).dots(channel).channels = [locDots(IDs,1:2) ones(length(locDots(IDs,1)),1)];
        end
    end
    toc
    
    
    
    %% divide into cells, then run analysis as normal
    ROI_path = [loc_pat filesep '..'];
    roi_prefix = 'RoiSet_Pos';
    fullpath = [ROI_path filesep roi_prefix num2str(my_pos) ];
    vertex = selfseg(fullpath);
    
    
    
    for cell_I = 1:size(vertex,2)
        
        %% keep only dots in THIS segmented cells
        fprintf('Selecting dots in this cell ...\n')
        for bcrI = 1:hyb
            for channel = 1:channum
                loc_dots = points(bcrI).dots(channel).channels(:,1:2);
                keep_I = [];
                keep_I = inpolygon(loc_dots(:,1),loc_dots(:,2),vertex(cell_I).x,vertex(cell_I).y);
                thisCellPoints(bcrI).dots(channel).channels = points(bcrI).dots(channel).channels(keep_I,:);
            end
        end
        
        fprintf('Finding barcodes...\n')
        [cell_I my_pos]
        tic
        [foundbarcodes] = gene_calling(thisCellPoints,radius);
        toc

        fprintf('Handle ambiguous codes...\n')
        [cell_I my_pos]
        tic
        [foundbarcodes, totdropped, copynumfinal, rawfound] = ...
            handle_ambiguous_barcodes(foundbarcodes, thisCellPoints, hyb, channum, ...
            barcodekey, alloweddiff,conthresh, ROI_path, [0 0],my_pos,roi_prefix)
        toc
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Find gene location ...\n')
        [cell_I my_pos]
        tic
        [dotlocations,copynumfinalrevised,PosList] = ...
            pointLocations(ROI_path,my_pos,hyb, channum, thisCellPoints, foundbarcodes,barcodekey,copynumfinal,radius)
        toc
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        fprintf('Counting seeds per gene ...\n')
        [cell_I my_pos]
        tic

        [seeds, int, ints] = numseeds(PosList,dotlocations);
        toc


        save([loc_pat filesep sprintf('%.0f_Analysis',channel_nm) filesep sprintf('Analysis_Details_NO_FISH_RCE_%.1f',TFperNom) filesep 'postProcData' filesep sprintf('Radius_%.1f_loc',radius) filesep  ...
            sprintf('Pos%.0f_Cell_%.0f_%.0fnm_results.mat',my_pos,cell_I,channel_nm)],'foundbarcodes','totdropped','copynumfinal',...
            'dotlocations','copynumfinalrevised','PosList',...
            'seeds','rawfound')
	    clear foundbarcodes rawfound
    end
end
