function [tforms] = alltforms(PathName,colors,channels,reference)
%   colors = # of channels in stack minus DAPI
%   channels = channels numbers that need to be aligned
%   reference = number of reference channel

fld = pwd;
Miji;
cd(fld);

if isempty(PathName) 
    [FileName,PathName,FilterIndex] = uigetfile('.tif');
    path_to_fish = ['path=[' PathName FileName ']'];
else
    path_to_fish = ['path=[' PathName ']'];
    C = strsplit(PathName,'\');
    FileName = C{end};
    k = strfind(PathName,'\');
    PathName = PathName(1:k(end)-1);
end

MIJ.run('Open...', path_to_fish);
%range = inputdlg('Where Should I Start In Z');
MIJ.run('Z Project...', 'projection=[Max Intensity]');
%MIJ.run('Subtract Background...', 'rolling=10 stack');
MIJ.run('Split Channels');
mkdir([PathName 'channelsWBG']);
for i = 1:colors
    name = ['C' num2str(i) '-MAX_' FileName];
    Images = uint16(MIJ.getImage(name));
    saveastiff(Images, [PathName 'channelsWBG\channel' num2str(i) '.tif'])
end

MIJ.run('Close All');
MIJ.exit;
for i = 1:colors
    tforms{i} = [];
end
for i = 1:length(channels)
    [tforms{channels(i)}] = findtformV3([PathName 'channelsWBG\channel' num2str(channels(i)) '.tif'],[PathName 'channelsWBG\channel' num2str(reference) '.tif']);
    pause('on');
    pause;
    close all;
end

tforms(cellfun(@isempty,tforms)) = tforms(reference);


