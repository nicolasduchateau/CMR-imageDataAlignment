clear all;
close all;
clc;

addpath('./Functions');

folderRoot = '../DATA/';

%% Read EXCEL file
ExcelFileName = [folderRoot,'0_ALL.xlsx'];
[num,txt,raw] = xlsread(ExcelFileName,1);

numCases =  size(raw,1)-1;

isValid = [[raw{2:end,12}]',[raw{2:end,13}]'];

%% Routine to identify (1) the LV junction + valve openings, (2) the basal and apical slices

selectedSlices_ALL = NaN(numCases,2,2);
sliceNumberIN_ALL = NaN(numCases,2);
for ki=1:numCases
    k = ki;
    disp(k);
    
    selectedSlices = NaN(2,2);
    for mi=1:2
        if isValid(ki,mi)==1
            selectedSlices(mi,:) = [raw{ki+1,6+(mi-1)*2},raw{ki+1,7+(mi-1)*2}];
            selectedSlices_ALL(ki,mi,:) = selectedSlices(mi,:);
        end
    end
    preprocessCVI42Data(k,folderRoot,selectedSlices);
    
    tmpmatD = [folderRoot , 'PROCESSED/D',sprintf('%03.0f',k),'.mat'];
    if exist(tmpmatD,'file')==2
        load(tmpmatD,'DCM');
        for mi=1:2
            if isValid(ki,mi)==1
                tmpL = length( DCM{mi}.USliceLocation );
                tmpDP = selectedSlices_ALL(ki,mi,2) - tmpL;
                if tmpDP < 0
                    tmpDP = 0;
                end
                tmpDM = 1 - selectedSlices_ALL(ki,mi,1);
                if tmpDM < 0
                    tmpDM = 0;
                end
                sliceNumberIN_ALL(ki,mi) = (tmpDP + tmpDM);
            end
        end
    end
    
end

