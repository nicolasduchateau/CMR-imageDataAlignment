function preprocessCVI42Data(k,folderRoot,selectedSlices)

tmpmatD = [folderRoot , 'PROCESSED/D',sprintf('%03.0f',k),'.mat'];

fileData.folderRoot = folderRoot;

if exist(tmpmatD,'file')==2

    load(tmpmatD,'CONTOURS');
    load(tmpmatD,'DCM');

    %%% Gets both MVO and LATE for each case
    [data_img,data_myo,data_inf,data_mvo]= cvi2Matrix(CONTOURS,DCM);

    originList = cell(2,1);
    myoOpeningList = cell(2,1);
    typeList = {'EGE','LGE'};

    for mi=1:2
        
        if sum( ~isnan(selectedSlices(mi,:)) ) > 0

            %%% Set the origin = the LV-RV junction
            fileData.fileName = [folderRoot,'P_DATA/ORIGIN/','Origin_',typeList{mi},'_',sprintf('%03.0f',k),'.mat'];
            if exist(fileData.fileName,'file') == 2
                load(fileData.fileName,'origine');
                load(fileData.fileName,'myoOpening');
                originList{mi} = origine;
                myoOpeningList{mi} = myoOpening;
            else
                if mi == 1 %%% no lesion drawn
                    infT = zeros(size(data_inf{mi}));
                else
                    infT = data_inf{mi};
                end
                [originList{mi},myoOpeningList{mi}] = initializeOrigin(data_img{mi},data_myo{mi},infT,data_mvo{mi},fileData,CONTOURS{mi},DCM{mi},selectedSlices(mi,:));
            end

            %%% Define the local coordinates + the segmented image
            fileData.fileName = ['CoordsAndData_',typeList{mi},'_',sprintf('%03.0f',k),'.mat'];
            fileData.fileNameMAT = [fileData.folderRoot,'P_DATA/COORDS/',fileData.fileName];

%             if exist(fileData.fileNameMAT,'file') ~= 2
            if (1)
                samplingFactor = 4;
                computeLocalCoords(originList{mi},myoOpeningList{mi},selectedSlices(mi,:),fileData,CONTOURS{mi},DCM{mi},samplingFactor);
            end
        end
        
    end
end

end
