function corrections = back_correct_multi_pos(blankIm)
%% Use the blank slide image to do corrections
cy5back = [];
bic = 0;
for posI = 1:length(blankIm)
    posI
    for i = 1:size(blankIm{posI},3)
        % blur blank ims
        bic = bic+1
        cy5back(:,:,bic) = imopen(blankIm{posI}(:,:,i),strel('disk',100));
    end
end
% find median over all blank im for each pixel
cy5med = median(cy5back,3);
% normalize to max
corrections = double(cy5med)/double(max(max(cy5med)));


