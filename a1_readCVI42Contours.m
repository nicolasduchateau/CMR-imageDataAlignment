clear all;
close all;
clc;

addpath('./FileTools');
addpath('./Functions');

folderRoot = '../DATA/';

%%% Read EXCEL file
ExcelFileName = [folderRoot,'0_ALL.xlsx'];
[num,txt,raw] = xlsread(ExcelFileName,1);

meshList = {'epi','endo'};
typeList = {'EGE','LGE'};

plotOUTPUT = 1;

numCases = size(num,1);
isValid = zeros(numCases,2);
for pi=1%1:numCases

    clear s;
    clear XML;
    clear DCM;
    
    tmpfile42 = [folderRoot , 'EXPORT/W',sprintf('%03.0f',pi),'.cvi42wsx'];
    tmpmat  = [folderRoot , 'PROCESSED/W',sprintf('%03.0f',pi),'.mat'];
    tmpmatD = [folderRoot , 'PROCESSED/D',sprintf('%03.0f',pi),'.mat'];

    if ( raw{pi+1,11}==1 ) %%% 0/1 analyze this case or not
        
        OKtoPlot = 0;

        %%% Read already extracted data
        if( exist(tmpmat,'file')==2 && exist(tmpmatD,'file')==2 )
            
            disp(['#',num2str(pi),' - Reading: ',tmpmat]);            
            load(tmpmat,'s');
            
            disp(['#',num2str(pi),' - Reading: ',tmpmatD]);
            load(tmpmatD,'CONTOURS');
            load(tmpmatD,'DCM');
            
            OKtoPlot = 1;
                
        else
        %%% Extract data from the DICOM
            
            if ( exist(tmpfile42,'file')==2 )

                disp(['#',num2str(pi),' - Reading DICOM data...']);

                DCM = cell(2,1);
                for mi=1:2 %%% 1=EGE / 2=LGE

                    %%% Parse DCMfolder and get UID
                    if (~isnan(raw{pi+1,2+mi}))
                        isValid(pi,mi) = 1;
                        str  = raw{pi+1,2};
                        if str(1) == '$'
                            str = strrep(str,'$','+');
                        end
                        DCM{mi}.folder = [folderRoot , 'DICOM/',sprintf('%03.0f',pi),'/',str,'/',raw{pi+1,2+mi},'/'];
                        listing = dir(DCM{mi}.folder);
                        listing([1,2]) = [];
                        numSlices = length(listing);
                        DCM{mi}.UID = cell(numSlices,1);
                        DCM{mi}.USliceLocation = zeros(numSlices,1);
                        for sci=1:numSlices
                            tmpfile = [DCM{mi}.folder,listing(sci).name];
                            tmpI = dicominfo(tmpfile);
                            if(sci == 1)
                                DCM{mi}.data = zeros(tmpI.Height,tmpI.Width,numSlices);
                                DCM{mi}.ImagePositionPatient = tmpI.ImagePositionPatient;
                                DCM{mi}.offset = zeros(1,3);
                                DCM{mi}.offset(3) = tmpI.SliceLocation;
                                DCM{mi}.spacing = tmpI.PixelSpacing;
                                DCM{mi}.spacing(3) = tmpI.SliceThickness;
                            end
                            DCM{mi}.data(:,:,sci) = dicomread(tmpfile);
                            DCM{mi}.UID{sci} = tmpI.MediaStorageSOPInstanceUID;
                            DCM{mi}.USliceLocation(sci) = tmpI.SliceLocation;
                        end

                        %%% check slices order from valve to apex
                        if mean( DCM{mi}.USliceLocation )>0 && raw{pi+1,5}==0 %%% 4+mi
                            disp('Inverted slices');
                            tmp_data = DCM{mi}.data;
                            for sci=1:numSlices
                                slicePos = numSlices - sci + 1;
                                tmp_data(:,:,sci) = DCM{mi}.data(:,:,slicePos);
                            end
                            DCM{mi}.data = tmp_data;
                            DCM{mi}.UID = flipud(DCM{mi}.UID);
                            DCM{mi}.USliceLocation = -flipud(DCM{mi}.USliceLocation);
                            DCM{mi}.spacing(3) = -DCM{mi}.spacing(3);
                        end

                    else
                        disp(['#',num2str(pi),' - Missing: ',typeList{mi}]);
                    end
                end

                disp(['#',num2str(pi),' - Converting XML to STRUCT...']);
                tic;
                s = xml2struct_S(tmpfile42);
                toc;

                disp(['#',num2str(pi),' - Writing: ',tmpmat]);
                save(tmpmat,'s');
                
                disp(['#',num2str(pi),' - Analyzing XML...']);
                
                XML = analyzeXML(s);

                disp(['#',num2str(pi),' - Extracting contours...']);

                CONTOURS = getContours(XML,DCM); %%% EGE/LGE - slice,contourType - nbContour - xyz coords

                %%% crop DCM to save space
                for mi=1:2
                    if(~isempty(CONTOURS{mi}))
                        tmpC = [];
                        for sci=1:size(CONTOURS{mi}.data,1)
                            for ci=1:2%size(CONTOURS{mi}.data,2) %%% EPI-ENDO
                                tmp = CONTOURS{mi}.data{sci,ci};
                                for r=1:length(tmp) %%% if several contours
                                    tmpC = [tmpC ; tmp{r}(:,1:2)];
                                end
                            end
                        end
                        if ~isempty(tmpC)
                            minT = min(tmpC); %%% both dims at once
                            maxT = max(tmpC); %%% both dims at once
                            minT = (minT - DCM{mi}.offset(1:2))./DCM{mi}.spacing(1:2)' + 1;
                            maxT = (maxT - DCM{mi}.offset(1:2))./DCM{mi}.spacing(1:2)' + 1;
                            minT = round(minT) - 10;
                            maxT = round(maxT) + 10;
                            DCM{mi}.data = DCM{mi}.data(minT(2):maxT(2),minT(1):maxT(1),:); %%% i,j

                            minT = (minT-1).*DCM{mi}.spacing(1:2)' + DCM{mi}.offset(1:2);
                            maxT = (maxT-1).*DCM{mi}.spacing(1:2)' + DCM{mi}.offset(1:2);
                            DCM{mi}.offset(1) = minT(1); %%% x
                            DCM{mi}.offset(2) = minT(2); %%% y
                        end
                    end
                end

                disp(['#',num2str(pi),' - Writing: ',tmpmatD]);
                save(tmpmatD,'CONTOURS');
                save(tmpmatD,'DCM','-append');
                
                OKtoPlot = 1;
                
            else        
                disp(['#',num2str(pi),' - Missing cvi42wsx file: ',tmpfile42]);
                
                OKtoPlot = 0;
                
            end
        end
    
        if OKtoPlot && plotOUTPUT
            plotContours(CONTOURS,DCM);
        end
    end
%     close all;
end
