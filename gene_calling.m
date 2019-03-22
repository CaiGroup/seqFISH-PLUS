function [foundbarcodes] = gene_calling(points, radius)

    hyb = 4;
    channum = 20;
    channels = channum;
    debug = 0;

    disp('Initializing....')
    %initialize barcode matrix
    for i = 1:hyb
        for j = 1:channum
            c = size(points(i).dots(j).channels,1);
            filler = zeros(c,1);
            filler(:) = j; 
            idxfill = zeros(c,1);
            idxfill(:) = 1:c;
            foundbarcodes(i).found(j).channel(:,i) = num2cell(filler);
            foundbarcodes(i).found(j).idx(:,i) = num2cell(idxfill);
        end
    end
    for i = 1:hyb
        for k = 1:hyb
            if k ~= i
                [i, k]
                for j = 1:channum
                        calledall = [];
                        calledidx = [];
                        for l = 1:channum
                                [~, ~, called1, ~, pair ] = colocalizeBarcodeV2_markedup_NO_INTENSITY(points, i,k, [j l],radius,debug );
                                calledall = [calledall, l*double(called1)];
                                brat = zeros(length(called1),1);
                                if ~isempty(pair)
                                    brat(pair(:,1)) = pair(:,2);
                                end
                                calledidx = [calledidx, brat];
                        end
                        calledall = mat2cell(calledall, ones(length(called1),1),channum);
                        x = cellfun(@removezeros,calledall,'UniformOutput',0);
                        foundbarcodes(i).found(j).channel(:,k) = x;
                        calledidx = mat2cell(calledidx, ones(length(called1),1),channum);
                        foundbarcodes(i).found(j).idx(:,k) = calledidx;
                end
            end
        end
    end
end
