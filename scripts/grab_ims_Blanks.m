function [My_Image] = grab_ims_Blanks(Path, posList, myChannel, zSlices)
warning ('off','all');
% load raw images
for posI = 1:length(posList)
    myPos = posList(posI)
    position = sprintf('Pos%.0f',myPos);

    My_Image(posI) = {loadTiffStacks(Path,position,(1:zSlices)+myChannel*zSlices)};
end
warning ('on','all');
