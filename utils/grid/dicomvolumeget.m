% Usage: 
 
function [header , volume]=dicomvolumeget(dirname)


find_file=1;

cd (dirname)
filelist=dir;% find out how many

[totalims,dum]=size(filelist);
totalims=totalims-2;
checkfile=zeros(totalims,1);
for target=1:totalims,
    find_file=1;
    target
    if find_file==1
        for test=0:totalims,
            if checkfile(test+1)==0
                filename2=filelist((test)+3,:).name;%1st two files are . and ..  (i.e. junk) so skip them
                
                imagefilenametest=  sprintf('%s%s',dirname,filename2);
                info=dicominfo(imagefilenametest);
                if info.InstanceNumber==target
                    imagefilename=imagefilenametest;
                    %imagefilename
                    checkfile(test+1)=1;
                    break;
                end
            else
                %do nothing
            end
        end
        find_file=0;
    end
    if target==1
        test=dicomread(imagefilename);
        [Nx,Ny]=size(test);
        volume=zeros(Nx,Ny,totalims);
    end
    imagex=dicomread(imagefilename);
    volume(:,:,target)=double(imagex);
end

header=dicominfo(imagefilename);











