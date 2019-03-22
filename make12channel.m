function offsets = make12channel(PathName, FolderName, hybs, channelstotal, channelsperhyb)

for i = 1:hybs
    listing{i} = dir([PathName '\' num2str(i) '\*.tif']);
end
mkdir(PathName,FolderName);
fld = pwd;
Miji;
cd(fld);
k2 = strfind(listing{1}(1).name, 'Pos');
        
for k = 1:length(listing{1})
    roi = listing{1}(k).name(k2(1):k2(1)+4);
    mkdir([PathName '\' FolderName],roi);
    for i = 0:(hybs/(channelstotal/channelsperhyb))-1
        counter = 1;
        for j = (i*channelstotal/channelsperhyb)+1:channelstotal/channelsperhyb*(i+1)
            Path = ['path=[' PathName '\' num2str(j) '\' listing{j}(k).name ']'];
            MIJ.run('Open...', Path);
            hybtemp = uint16(MIJ.getCurrentImage);
            if mod(j,channelstotal/channelsperhyb) ~=1
                regi = cat(3,regi,hybtemp(:,:,4));
%                 c = normxcorr2(regi(:,:,counter),regi(:,:,1));
%                 [~, imax] = max(abs(c(:)));
%                 [ypeak, xpeak] = ind2sub(size(c),imax(1));
%                 offset{i+1,counter} = [(xpeak-size(regi(:,:,1),2)) (ypeak-size(regi(:,:,1),1))];
%                 hybtempreg = imtranslate(hybtemp,offset{i+1,counter});
                tform = imregcorr(regi(200:800,200:800,counter),regi(200:800,200:800,1));
                hybtempreg = imwarp(hybtemp,tform,'OutputView',imref2d(size(regi(:,:,1))));
                hyb = cat(3,hyb, hybtempreg(:,:,1:3));
                regireg = cat(3,regireg, hybtempreg(:,:,4));
                offsets{i+1,counter} = tform;
            else
                hyb = hybtemp(:,:,1:3);
                regireg = hybtemp(:,:,4);
                regi = hybtemp(:,:,4);
            end
            counter = counter + 1;
            MIJ.run('Close All')
        end
        MIJ.createImage('HYB',hyb,true);
        MIJ.run('Stack to Hyperstack...', ['order=xyczt(default) channels=' num2str(channelstotal) ' slices=1 frames=1 display=Grayscale']);
        MIJ.run('Save', ['save=[' PathName '\' FolderName '\' roi '\' num2str(i+1) '.tif' ']']);
        MIJ.run('Close All')
        MIJ.createImage('Regireg',regireg,true);
        MIJ.run('Save', ['save=[' PathName '\' FolderName '\' roi '\hyb' num2str(i+1) 'RegistrationCheck.tif' ']']);
        MIJ.run('Close All')
    end
end
        
MIJ.exit        
        
        
        
        
        
        
        
  