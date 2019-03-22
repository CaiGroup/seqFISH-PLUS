function [foundbarcodes, totdropped, copynumfinal, rawfound] = ...
    handle_ambiguous_barcodes(foundbarcodes, points, hyb, channum, barcodekey, alloweddiff,conthresh, PathName, regvec, posnum, roi_prefix)

channels = 1:channum;

% adapted for Sheel's BarcodeNoMiji_v5.m
rawfound = foundbarcodes;
disp('Determining Barcodes....')

%call Barcodes
totdropped = 0;
% 1: per dot, due to not enough dots close enough in at least 2 other
%       rounds.
% 2: per dot, dropped because multi and not sufficiently close to the code
%       book
for k = 1:hyb
    for j = 1:channum
        drop = [];
        bobo = [];
        bratat = {};
        br = [];
        len = [];
        multi = [];
        rows = [];
        
        %lentemp: number of barcodes with this hyb and color as seed
        lentemp = size(foundbarcodes(k).found(j).channel,1);
        for i = 1:lentemp
            bobo = cell2mat(foundbarcodes(k).found(j).channel(i,:));
            % bratat contains all hyb possible channels concatinated to 1
            % array
            bratat{i,1} = bobo; 
        end
        % br: number of hybs with zero for one dot
        br = cellfun(@(x) sum(x == 0),bratat,'UniformOutput',0);
        %initialize called barcodes as the size of the dots
        called = zeros(lentemp,1);
        % more than 1 non called hyb...
        drop=cellfun(@(x) x>1,br);
            % set whole thing to zero, presumably to drop later
        foundbarcodes(k).found(j).channel(drop,:) ={0};
        foundbarcodes(k).found(j).idx(drop,:) ={0};
        bratat(drop) = {[0 0 0 0]};
        br2 = cellfun(@removezeros2,bratat,'UniformOutput',0);
        % br2 is bratat without zeros
        % len is the non-zero length of each baratat entry
        len = cellfun(@length,br2);
        multi = len+cell2mat(br)>hyb;
        [rows] = find(multi ==1);
        % deal with mulitple hit codes
        for l = 1:length(rows)
            possbarcodes = [];
            Dposs = [];
            A = [];
            posscell = [];
            code = [];
            gene = [];
            possreal = [];
            vernoi = [];
            vernoiz = [];
            ind = [];
            set = [];
            vari = [];
            variz = [];
            % find all equally minimum hamming distance to the code book
            % combinations of this dot
            possbarcodes = combvec(foundbarcodes(k).found(j).channel{rows(l),:})';
            Dposs = pdist2(possbarcodes,barcodekey.barcode,'hamming');
            A = Dposs == min(min(Dposs));
            % code is the idea of the posscodes that are closest
            [code, ~] = find(A == 1);
            % if there is only one possible code and it is close enough
            % call that the correct one
            if length(code) == 1 && Dposs(A)*hyb < alloweddiff
                posscell = num2cell(possbarcodes);
                % set found barcode to the closest barcode
                foundbarcodes(k).found(j).channel(rows(l),:) = posscell(code,:);
                bratat{rows(l)} = cell2mat(posscell(code,:));
            elseif Dposs(A)*hyb < alloweddiff & length(code) > 1
                % grab all the equally closest ones
                possreal = possbarcodes(code,:);
                % idx contains the index for each channel in each hyb that
                % showed up for the dot seeded from channel j in hyb k
                ind = foundbarcodes(k).found(j).idx(rows(l),:);
                caca = zeros(1,length(channels));
                caca(1,j) = ind{k};
                ind{k} = caca;
                for p = 1:size(possreal,1)
                    for h = 1:hyb
                        % for the non-zero values in a given possible real
                        % code
                        if possreal(p,h) > 0
                            % for each possible code, build a location list
                            % of the dots associated with it
                            set(p).set(h,:) = points(h).dots(possreal(p,h)).channels(ind{h}(possreal(p,h)),:);
                        else 
                            set(p).set(h,:) = zeros(1,3);
                        end
                    end
                    % all the dots should have the same location
                    % so find the sum of the variance in position
                    vari(p) = sum(var(set(p).set));
                    variz(p) = var(set(p).set(:,3));
                end
                % binary map of the minimum variance codes
                vernoi = vari == min(vari);
                vernoiz = variz == min(variz);
                if sum(vernoi) == 1
                    foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(vernoi',:));
                    bratat{rows(l)} = possreal(vernoi',:);
                end
            else
                foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                bratat(rows) = {[0 0 0 0]};
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fill in dropped cells
        % br2 is bratat (cell array containg barcode arrays (as opposed
        % cellarry barcodes)) with zeros removed from the barcode arrays
        br2 = cellfun(@removezeros2,bratat,'UniformOutput',0);
        len = cellfun(@length,br2);
        % find the codes tha are one hyb value less than the full number
        missing = len == hyb-1;
        % rows contains their IDs
        [rows] = find(missing ==1);
        for l = 1:length(rows)
            possbarcodes = [];
            Dposs = [];
            A = [];
            posscell = [];
            code = [];
            gene = [];
            possreal = [];
            vernoi = [];
            vernoiz = [];
            ind = [];
            set = [];
            vari = [];
            variz = [];
            % essentially a repeat of above, but this time fill in the
            % missing value to make a four round code
            possbarcodes = combvec(foundbarcodes(k).found(j).channel{rows(l),:})';
            Dposs = pdist2(possbarcodes,barcodekey.barcode,'hamming');
            A = Dposs == min(min(Dposs));
            [code, gene] = find(A == 1);
            if length(gene) == 1 && Dposs(A)*hyb < alloweddiff
                foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(barcodekey.barcode(gene,:));
                bratat{rows(l)} = barcodekey.barcode(gene,:);
            % if there is more than one that are equally close and close
            % enough
            elseif Dposs(A)*hyb < alloweddiff & length(code) > 1
                possreal = possbarcodes(code,:);
                ind = foundbarcodes(k).found(j).idx(rows(l),:);
                caca = zeros(1,5);
                caca(1,j) = ind{k};
                ind{k} = caca;
                for p = 1:size(possreal,1)
                    for h = 1:hyb
                        if possreal(p,h) > 0
                            set(p).set(h,:) = points(h).dots(possreal(p,h)).channels(ind{h}(possreal(p,h)),:);
                        else 
                            set(p).set(h,:) = zeros(1,3);
                        end
                    end
                    vari(p) = sum(var(set(p).set));
                    variz(p) = var(set(p).set(:,3));
                end
                vernoi = vari == min(vari);
                vernoiz = variz == min(variz);
                % grab the one with smallest spatial variance in dots
                % called
                if sum(vernoi) == 1
                    foundbarcodes(k).found(j).channel(rows(l),:) = num2cell(possreal(vernoi',:));
                    bratat{rows(l)} = possreal(vernoi',:);
                else
                        totdropped = totdropped + 1;
                        foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                        foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                        bratat(rows) = {[0 0 0 0]};
                end
            else
                foundbarcodes(k).found(j).channel(rows(l),:) ={0}; % remove ambiguous dots
                foundbarcodes(k).found(j).idx(rows(l),:) ={0};
                bratat(rows(l)) = {[0 0 0 0]};
            end
        end 
        
        
        % call barcodes
        if ~isempty(foundbarcodes(k).found(j).channel)
            cell2 = [];
            minmat = [];
            tester = [];
            logicalstuff = [];
            posscell = {};
            % initiallizes posscell to add barcodes seeded from this 
            % hyb (bcR, k) and channel (hyb, j)
            [posscell{1:size(foundbarcodes(k).found(j).channel,1),1:hyb}] = deal(0);
            % finds distance from each identified code to each code 
            D = pdist2(cell2mat(foundbarcodes(k).found(j).channel),barcodekey.barcode,'hamming');
            [r,c] = size(D);
            minmat = cellfun(@(x) x == min(x,[],2) & x<=((alloweddiff-1)/hyb),mat2cell(D,ones(1,r),c),'UniformOutput',0);
            % finds the ones from the previous step, i.e. barcodes to call
            [re, co] = cellfun(@(x) find(x),minmat,'UniformOutput',0);
            dasdf = cellfun(@isempty,re,'UniformOutput',0);
            re(logical(cell2mat(dasdf))) = {0};
            mm = cellfun(@length,re);
            mmdrop = mm > 1;
            re(mmdrop) = {0};
            % finds the cells with a code (i.e. re == 1)
            rows = find(cell2mat(re));
            % populates posscell with found code book codes
            posscell(rows,:) = num2cell(barcodekey.barcode(cell2mat(co(rows)),:));
            co(logical(cell2mat(dasdf))) = {0};
            co(mmdrop) = {0};
            foundbarcodes(k).found(j).called = cell2mat(co);
            awe = foundbarcodes(k).found(j).called > 0;
            rawfound(k).found(j).called = foundbarcodes(k).found(j).called;
            for m = 1:size(foundbarcodes(k).found(j).channel,1)
                % for each hyb (barcode number)
                for n = 1:hyb
                    mu = zeros(1,channum);
                    if foundbarcodes(k).found(j).channel{m,n} > 0
                        mu(1,foundbarcodes(k).found(j).channel{m,n}) = 1;
                    end
                    cell2{m,n} = mu;
                end
            end
            foundbarcodes(k).found(j).idx = cellfun(@(x,y) x.*y,foundbarcodes(k).found(j).idx,cell2,'UniformOutput',0);
        else
            foundbarcodes(k).found(j).idx = [];
        end
    end
end


disp('Forming Point Consensus....')
%consensus point calling
for i = 1:hyb
    for j = 1:length(channels)
        % initialize compiled with zeros the size of the number of seeded
        foundbarcodes(i).found(j).compiled = zeros(size(foundbarcodes(i).found(j).called,1),hyb);
    end
end
% create a compiled list for every dot
disp('Compiling List....')
% compiled will contain an with the code called using that dot based on
% seeding from each hyb, which will be used to form consensus of where each
% dot belongs
for i = 1:hyb
    for j = 1:length(channels)
        calledrows = find(foundbarcodes(i).found(j).called>0);
        if ~isempty(calledrows)
            % r contains cells array, seeded x hybs, with ones for kepth
            % codes
            % c contains the a cell array with the channel for each hyb for
            % each kept code
            % v contains the a cell array with the point id associated with
            % each value in c
            [r,c,v]=cellfun(@(x) find(x),foundbarcodes(i).found(j).idx,'UniformOutput',0);
            for k = 1:length(calledrows)
                for l = 1:hyb
                    if ~isempty(c{calledrows(k),l})
                        foundbarcodes(l).found(c{calledrows(k),l}).compiled(v{calledrows(k),l},i) = foundbarcodes(i).found(j).called(calledrows(k));
                    end
                end
            end
        end
    end
end


disp('Dropping Ambiguous Matches....')
for i = 1:hyb
    for j = 1:length(channels)
        % Drop Ambigous Matches
        % conver to compiled to cell array
        cellver = mat2cell(foundbarcodes(i).found(j).compiled,ones(size(foundbarcodes(i).found(j).compiled,1),1),hyb);
        % remove zeros, so they can not be the mode
        cellver = cellfun(@removezeros, cellver,'UniformOutput',0);
        % resolve ambiquous by grabbing the most frequently called code
        [M,F,C] = cellfun(@mode,cellver,'UniformOutput',0);
        % if there are multiple mode values M only contains the smalles,
        % whereas C contains all of them ...
        eq = cellfun(@(x) length(x{1}),C);
        M = cell2mat(M);
        % if there's more than one mode, then toss it
        M(eq>1)=0;
        % F says how often that code was called, so this makes sure it was
        % called seeding from a sufficient number of hybs (i.e. >conthresh)
        M(cell2mat(F)<conthresh) = 0;
        % save the consensus code
        foundbarcodes(i).found(j).consensus = M;
    end
end


disp('Assigning to Cell....')
	[PathName filesep roi_prefix num2str(posnum) ]
    fullpath = [PathName filesep roi_prefix num2str(posnum) ];
    vertex = selfseg(fullpath);
    copynumfinal(:,:) = barcodekey.names;
    for i = 1:length(vertex)
        for j = 1:hyb
            %j = 1;
            allcalled = [];
            for k = 1:channum
                include(j).points(k).channel = inpolygon(points(j).dots(k).channels(:,1),points(j).dots(k).channels(:,2),vertex(i).x+regvec(1),vertex(i).y+regvec(2));
                allcalled = [allcalled; foundbarcodes(j).found(k).consensus(include(j).points(k).channel,:)];
            end
            copy = histc(allcalled(:),0:length(barcodekey.names));
            copynum(:,j) = copy(2:end);
        end
        copynumfinal(:,i+1) = num2cell(max(copynum,[],2));
        copynumfinalsum(:,i+1) = num2cell(sum(copynum,2));
    end

