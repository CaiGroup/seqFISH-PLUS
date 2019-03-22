function [images] = loadTiffStacks(path,position,channels)
% find the name of the tif stack
allFiles = dir(path);
for i = 1:length(allFiles)
    if ~isempty(strfind(allFiles(i).name,position)) && ~isempty(strfind(allFiles(i).name,'ome'))
        myFile = allFiles(i).name;
        break
    end
end
% open tif image for this position
FileTif = [path filesep myFile];
InfoImage = imfinfo(FileTif);
NumberImages = length(InfoImage);
TifLink = Tiff(FileTif,'r');


for i = 1:length(channels)
    channel = channels(i);
    TifLink.setDirectory(channel);
    images(:,:,i) = TifLink.read();
end
