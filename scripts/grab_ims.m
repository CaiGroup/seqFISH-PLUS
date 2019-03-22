function [My_Image] = grab_ims(Path, posList, numRnds, myChannel, zSlices)

warning ('off','all');
% load raw images
for posI = 1:length(posList)
    myPos = posList(posI);
    position = sprintf('Pos%.0f',myPos);
    % grab all rounds, bead only is first
    for i = 1:numRnds
        myPath = [Path filesep '..' filesep sprintf('HybCycle_%.0f',i-1)];
        My_Image(i) = {loadTiffStacks(myPath,position,(1:zSlices)+myChannel*zSlices)};
    end
end

