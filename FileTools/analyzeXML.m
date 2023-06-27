function OUT = analyzeXML(s)
% % source: Nicolas DUCHATEAU, CREATIS - Université Lyon 1

for hi=1:size(s.Workspaceu_colonu_Workspace.Hashu_colonu_item,2)
    if isfield(s.Workspaceu_colonu_Workspace.Hashu_colonu_item{1,hi},'Listu_colonu_item')
        tmpH = s.Workspaceu_colonu_Workspace.Hashu_colonu_item{1,hi}.Listu_colonu_item{1,1}.Hashu_colonu_item;
        break;
    end
end

for hi=1:size(tmpH,2)
    if isfield(tmpH{1,hi},'Hashu_colonu_item')
        tmpS = tmpH{1,hi};
        break;
    end
end

numStates = str2double(tmpS.Attributes.Hash_colon_count);
imageStates = cell(numStates,1);
for si=1:numStates  
    imageStates{si}.UID = tmpS.Hashu_colonu_item{1,si}.Attributes.Hash_colon_key;
    
    tmpH = tmpS.Hashu_colonu_item{1,si}.Hashu_colonu_item;
    for hi=1:size(tmpH,2)
        if strcmp( tmpH{1,hi}.Attributes.Hash_colon_key , 'Contours' )
            tmpC = tmpH{1,hi};
        end
    end
    
    numContours = str2double(tmpC.Attributes.Hash_colon_count);
    contours = cell(numContours,1);

    for ci=1:numContours
        if(numContours == 1)
            tmpCC = tmpC.Hashu_colonu_item;
        else
            tmpCC = tmpC.Hashu_colonu_item{1,ci};
        end
        contours{ci}.name = tmpCC.Attributes.Hash_colon_key;
        
        tmpH = tmpCC.Hashu_colonu_item;
        for hi=1:size(tmpH,2)
            if strcmp( tmpH{1,hi}.Attributes.Hash_colon_key , 'SubpixelResolution' )
                SbResolution = str2double(tmpH{1,hi}.Text);
                break;
            end
        end
        
        for hi=1:size(tmpH,2)
            if strcmp( tmpH{1,hi}.Attributes.Hash_colon_key , 'Points' )
                tmpHP = tmpCC.Hashu_colonu_item{1,hi};
                break;
            end
        end
        
        numPoints = str2double(tmpHP.Attributes.List_colon_count);
        points = zeros(numPoints,2);
        for pi=1:numPoints
            if numPoints == 1
                tmpP = tmpHP.Listu_colonu_item;
            else
                tmpP = tmpHP.Listu_colonu_item{pi};
            end
            points(pi,1) = str2double(tmpP.Pointu_colonu_x.Text);
            points(pi,2) = str2double(tmpP.Pointu_colonu_y.Text);
        end
        contours{ci}.points = points / SbResolution;
    end
    imageStates{si}.contours = contours;

end
OUT.imageStates = imageStates;

end
