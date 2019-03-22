

PathName = uigetdir;

listing = dir([PathName '/*.tif']);
fld = pwd;
Miji;
cd(fld);    
cy7back = [];
cy5back = [];
a594back = [];
cy3back = [];
a488back = [];

for i = 1:length(listing)
    path_to_fish = ['path=[' PathName '\' listing(i).name ']'];
    MIJ.run('Open...', path_to_fish);
    MIJ.run('Split Channels');
    for loop = 1:3
        name = ['C' num2str(loop) '-' listing(i).name];
        im = uint16(MIJ.getImage(name));
        hybnum.color{loop} = im;
    end
    MIJ.run('Close All');
    %imcortex = loadtiff(['/Volumes/New Volume/10132015Brain/Final Images to Use/Organized/pos' num2str(num(i)) '/Nissel.tif']);
    %S = load([PathName '/pos' num2str(num(i)) '/hybnum' num2str(num(i)) '.mat']);
    cy7 = max(hybnum.color{1},[],3);
    cy5 = max(hybnum.color{2},[],3);
    a594 = max(hybnum.color{3},[],3);
%     cy3 = max(hybnum.color{4},[],3);
%     a488 = max(hybnum.color{5},[],3);
    cy7back(:,:,i) = imopen(cy7,strel('disk',100));
    cy5back(:,:,i) = imopen(cy5,strel('disk',100));
    a594back(:,:,i) = imopen(a594,strel('disk',100));
%     cy3back(:,:,i) = imopen(cy3,strel('disk',100));
%     a488back(:,:,i) = imopen(a488,strel('disk',100));
end

cy7med = median(cy7back,3);
cy5med = median(cy5back,3);
a594med = median(a594back,3);
% cy3med = median(cy3back,3);
% a488med = median(a488back,3);

cy7med = double(cy7med)/double(max(max(cy7med)));
cy5med = double(cy5med)/double(max(max(cy5med)));
a594med = double(a594med)/double(max(max(a594med)));
% cy3med = double(cy3med)/double(max(max(cy3med)));
% a488med = double(a488med)/double(max(max(a488med)));

corrections(1).corrections{1} = cy7med;
corrections(1).corrections{2} = cy5med;
corrections(1).corrections{3} = a594med;
% corrections(1).corrections{4} = cy3med;
% corrections(1).corrections{5} = a488med;

figure;
surf(double(cy7med(1:8:end,1:8:end))),zlim([0 1]);
ax = gca;
ax.YDir = 'reverse';

figure;
surf(double(cy5med(1:8:end,1:8:end))),zlim([0 1]);
ax = gca;
ax.YDir = 'reverse';

figure;
surf(double(a594med(1:8:end,1:8:end))),zlim([0 1]);
ax = gca;
ax.YDir = 'reverse';

% figure;
% surf(double(cy3med(1:8:end,1:8:end))),zlim([0 1]);
% ax = gca;
% ax.YDir = 'reverse';
% 
% figure;
% surf(double(a488med(1:8:end,1:8:end))),zlim([0 1]);
% ax = gca;
% ax.YDir = 'reverse';


MIJ.exit
