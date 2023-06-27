function CONTOURS = getContours(XML,DCM)
% % source: Nicolas DUCHATEAU, CREATIS - Université Lyon 1

colorLabels = {'saepicardialContour','saendocardialContour','saReferenceMyoContour','saEnhancementReferenceMyoContour','noReflowAreaContour','excludeEnhancementAreaContour'};

CONTOURS = cell(2,1);
for mi=1:2 %%% MVO / LATE
    if(~isempty(DCM{mi}))
        CONTOURS{mi}.data = cell(size(DCM{mi}.data,3),length(colorLabels));
    end
end

for si=1:size(XML.imageStates,1)

    %%% find corresponding DCM
    for mi=1:2 %%% MVO / LATE
                
        colorFound = 0;
        sliceFound = 0;
        clear miOK sciOK colorID;
        if(~isempty(CONTOURS{mi}))
            for sci=1:length(DCM{mi}.UID)
                if( strcmp(XML.imageStates{si}.UID , DCM{mi}.UID{sci}) )
                    miOK = mi;
                    sciOK = sci;
                    sliceFound = 1;
                    break;
                end
                if(sliceFound == 1)
                    break;
                end
            end
        end
    
        if(sliceFound)
            for ci=1:size(XML.imageStates{si}.contours,1)
                tmp = XML.imageStates{si}.contours{ci};

                for cli=1:length(colorLabels)
                    if( strcmp(tmp.name,colorLabels{cli}) || strcmp(tmp.name(1:(end-4)),colorLabels{cli}) )
                        colorID = cli;
                        colorFound = 1;
                        break;
                    end
                end
                if(colorFound)
                    numPoints = size(tmp.points,1);
                    tmp.points(:,3) = repmat(DCM{mi}.USliceLocation(sciOK),[numPoints,1]);

                    tmpS = [DCM{miOK}.spacing(1),DCM{miOK}.spacing(2),1];
                    tmpH = DCM{miOK}.offset; tmpH(3) = 0;
                    numSubContours = size(CONTOURS{miOK}.data{sciOK,colorID},1);
                    CONTOURS{miOK}.data{sciOK,colorID}{numSubContours+1} = tmp.points .* repmat(tmpS,[numPoints,1]) + repmat(tmpH,[numPoints,1]);
                end
            end
        end
    end

end

end


