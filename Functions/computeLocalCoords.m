function computeLocalCoords(originListi,myoOpeningListi,selectedSlices,fileData,CONTOURSi,DCMi,samplingFactor)

[Nr,Nc,Nz] = size(DCMi.data);

%% READING DATA

colorList = {[0,1,0],[1,0,0],[0,0,1],[1,0,1],[0,1,1],[0.5,0.5,0.5]};
%%% 1=epi / 2=endo / 3=ROImyo / 4=ROIlesion / 5=lesion / 6=?

coords = NaN(Nr*samplingFactor,Nc*samplingFactor,Nz,3);
jetResolution = 512;
angleOUT = NaN(Nz,2); %%% min,max
coordsList = NaN(Nz,1);

data_IMG = NaN(Nr*samplingFactor,Nc*samplingFactor,Nz);
data_SEG = NaN(Nr*samplingFactor,Nc*samplingFactor,Nz);

for slice = 1:Nz
    sliceToDraw = slice; %%%%%%%%%%%%%%%%
    
    slicePos = slice;
    disp(['Slice ',num2str(slicePos),'/',num2str(Nz)]);

    if (~isempty(CONTOURSi.data{slicePos,1})) && (~isequal(originListi(slicePos,:),[0,0]))

        tmpI = squeeze( DCMi.data(:,:,slicePos) );
        tmpX = ([1,Nc]-1)*DCMi.spacing(1) + DCMi.offset(1);
        tmpY = ([1,Nr]-1)*DCMi.spacing(2) + DCMi.offset(2);
                        
        EPI  = CONTOURSi.data{slicePos,1}{1}(:,1:2);
        if ~isempty(CONTOURSi.data{slicePos,2})
            ENDOok = 1;
        else
            ENDOok = 0;
        end
            
        if ENDOok
            ENDO = CONTOURSi.data{slicePos,2}{1}(:,1:2);
        else
            ENDO = EPI;
        end

        datatmp = CONTOURSi.data{slicePos,1}; %%% EPI
        x1 = (datatmp{1}(:,1) - DCMi.offset(1))/DCMi.spacing(1)*samplingFactor + 1;
        y1 = (datatmp{1}(:,2) - DCMi.offset(2))/DCMi.spacing(2)*samplingFactor + 1;
        masktmp = poly2mask(x1,y1,Nr*samplingFactor,Nc*samplingFactor);
        maskEPI = masktmp;

        if ENDOok
            datatmp = CONTOURSi.data{slicePos,2}; %%% ENDO
            x1 = (datatmp{1}(:,1) - DCMi.offset(1))/DCMi.spacing(1)*samplingFactor + 1;
            y1 = (datatmp{1}(:,2) - DCMi.offset(2))/DCMi.spacing(2)*samplingFactor + 1;
            masktmp = poly2mask(x1,y1,Nr*samplingFactor,Nc*samplingFactor);
            maskENDO = masktmp;
        else
            maskENDO = maskEPI;
        end
        [I,J] = ind2sub([Nr*samplingFactor,Nc*samplingFactor],find(maskENDO==1));
        X = (J - 1)*DCMi.spacing(1)/samplingFactor + DCMi.offset(1);
        Y = (I - 1)*DCMi.spacing(2)/samplingFactor + DCMi.offset(2);
        centerENDO = [mean(X),mean(Y)]; %%% xy
        
        %%% remove 1 pixels alone
        se = strel('square',2*samplingFactor);
        maskEPI = imopen(maskEPI,se);
        maskENDO = imopen(maskENDO,se);
        if ENDOok
            maskMYO = maskEPI .* (~maskENDO);
        else
            maskMYO = maskEPI;
        end

        %%% upsampled image
        F = griddedInterpolant(tmpI);
        xq = linspace(1,Nr,Nr*samplingFactor);
        yq = linspace(1,Nc,Nc*samplingFactor);
        vq = F({xq,yq});
        vq(maskMYO==0) = NaN;
        data_IMG(:,:,slice) = vq;

        %%% upsampled segmentation / works directly on contours - no interp
        tmp = NaN(Nr*samplingFactor,Nc*samplingFactor);
        tmp(maskMYO==1) = 1;
        
        datatmp = CONTOURSi.data{slicePos,4}; %%% INFARCT
        for sbi=1:length(datatmp) %%% if several contours for a given structure (ex: MVO)
            if(~isempty(datatmp{sbi}))
                x1 = (datatmp{sbi}(:,1) - DCMi.offset(1))/DCMi.spacing(1)*samplingFactor + 1;
                y1 = (datatmp{sbi}(:,2) - DCMi.offset(2))/DCMi.spacing(2)*samplingFactor + 1;
                masktmp = poly2mask(x1,y1,Nr*samplingFactor,Nc*samplingFactor);
                tmp(masktmp==1) = 2;
            end
        end
                
        datatmp = CONTOURSi.data{slicePos,5}; %%% MVO
        for sbi=1:length(datatmp) %%% if several contours for a given structure (ex: MVO)
            if(~isempty(datatmp{sbi}))
                x1 = (datatmp{sbi}(:,1) - DCMi.offset(1))/DCMi.spacing(1)*samplingFactor + 1;
                y1 = (datatmp{sbi}(:,2) - DCMi.offset(2))/DCMi.spacing(2)*samplingFactor + 1;
                masktmp = poly2mask(x1,y1,Nr*samplingFactor,Nc*samplingFactor);
                tmp(masktmp==1) = 3;
            end
        end
        data_SEG(:,:,slice) = tmp;
        
        %%% contours and origin
        contourEPI  = double(bwperim(maskEPI));
        contourENDO  = double(bwperim(maskENDO));
        contourMYO = double(bwperim(maskMYO));
                
        tmpO = originListi(slicePos,:); %%% ij indices
        origine = (fliplr(tmpO) - 1).*DCMi.spacing(1:2)' + DCMi.offset(1:2); %%% xy

        %% define "forbidden" zone
        p1 = centerENDO; %%% xy
        p20 = squeeze( myoOpeningListi(slicePos,1,:) )';
        p20 = (fliplr(p20) - 1).*DCMi.spacing(1:2)' + DCMi.offset(1:2); %%% xy
        v20 = p20 - p1;
        p20new = p1 + 4*v20;
        p21 = squeeze( myoOpeningListi(slicePos,2,:) )';
        p21 = (fliplr(p21) - 1).*DCMi.spacing(1:2)' + DCMi.offset(1:2); %%% xy
        v21 = p21 - p1;
        p21new  = p1 + 4*v21;
        p22 = mean([p20;p21]);
        v22 = p22 - p1;
        p22new  = p1 + 10*v22;
        OPENold = [p1;p20;p22;p21];
        OPEN = [p1;p20new;p22new;p21new];
        
        if(slice == sliceToDraw)
            figure('units','normalized','position',[0,0,1,1],'Visible','off');
                    
            h1 = subplot(2,2,1); hold on;
            colormap(h1,gray);
            imagesc(tmpI,'XData',tmpX,'YData',tmpY);
            axis equal; grid on;
            axis ij;
            zoom(1); %%% zoom(3) if not cropped
            axis manual;

            plot(EPI(:,1),EPI(:,2),'g','LineWidth',2);
            plot(ENDO(:,1),ENDO(:,2),'r','LineWidth',2);
            for d=4:5
                tmp = CONTOURSi.data{slicePos,d};
                for d2=1:length(tmp)
                    tmp2 = CONTOURSi.data{slicePos,d}{d2}(:,1:2);
                    plot(tmp2(:,1),tmp2(:,2),'Color',colorList{d},'LineWidth',2);
                end
            end
            plot(origine(1),origine(2),'y.','MarkerSize',20);
            for d=[2,4]
                plot(OPENold(d,1),OPENold(d,2),'c.','MarkerSize',20);
                plot([OPEN(1,1),OPEN(d,1)],[OPEN(1,2),OPEN(d,2)],'c');
            end
            plot(centerENDO(1),centerENDO(2),'y+');

            title(['Slice ',num2str(slice)]);

            drawnow;            
        end
                
        %% RADIAL       
        
        p1 = centerENDO; %%% xy

        Inew = NaN(Nr*samplingFactor,Nc*samplingFactor);
        idx = find(maskMYO==1);
        for ti=1:length(idx)            
            %%% coordinates of each point of contour
            [I,J] = ind2sub([Nr*samplingFactor,Nc*samplingFactor],idx(ti));
            X = (J - 1)*DCMi.spacing(1)/samplingFactor + DCMi.offset(1);
            Y = (I - 1)*DCMi.spacing(2)/samplingFactor + DCMi.offset(2);
            p2 = [X,Y];
            
            if ~inpolygon(p2(1),p2(2),OPEN(:,1),OPEN(:,2))
                tmp = sum( (EPI - repmat(p2,[size(EPI,1),1])).^2 , 2 ); %%% L2 distance
                dToEPI = ( tmp( find(tmp==min(tmp),1) ) ).^(1/2);

                if ENDOok
                    tmp = sum( (ENDO - repmat(p2,[size(ENDO,1),1])).^2 , 2 ); %%% L2 distance
                    dToENDO = ( tmp( find(tmp==min(tmp),1) ) ).^(1/2);
                else
                    tmp = sum( (p1 - p2).^2 , 2 ); %%% L2 distance
                    dToENDO = ( tmp ).^(1/2);
                end

                Inew(I,J) = dToENDO / (dToENDO + dToEPI);
            end
        end
        coords(:,:,slice,1) = Inew;

        %% CIRCUMFERENTIAL
        
        p1 = centerENDO; %%% xy
        p20 = origine;
        v20 = p20 - p1;

        I1 = NaN(Nr*samplingFactor,Nc*samplingFactor);
        Inew = NaN(Nr*samplingFactor,Nc*samplingFactor);
        idx = find(maskMYO==1);
        for ti=1:length(idx)            
            %%% coordinates of each point of contour
            [I,J] = ind2sub([Nr*samplingFactor,Nc*samplingFactor],idx(ti));
            X = (J - 1)*DCMi.spacing(1)/samplingFactor + DCMi.offset(1);
            Y = (I - 1)*DCMi.spacing(2)/samplingFactor + DCMi.offset(2);
            p2 = [X,Y];            

            if ~inpolygon(p2(1),p2(2),OPEN(:,1),OPEN(:,2))
                v2 = p2-p1;
    %             angleT = atan2(v20(1),v20(2)) - atan2(v2(1),v2(2)); %%% (x,y) 
                angleT = atan2(v20(2),v20(1)) - atan2(v2(2),v2(1)); %%% (y,x)
                Inew(I,J) = angleT + 2*pi*(angleT<0);
            end
        end
        Inew = Inew ./ (2*pi);
        coords(:,:,slice,2) = Inew;
        
        %%% find missing angles
        t = linspace(0,1,360);
        tOUT = [];
        for ti=1:length(t)
            tmp = abs( Inew(:) - t(ti) );
            if(min(tmp) > 0.01)
                tOUT = [tOUT , t(ti)];
            end
        end
        if ~isempty(tOUT)
            if( tOUT(1)==0 && tOUT(end)==1 )
                tmp = tOUT(2:end) - tOUT(1:end-1);
                itmp = find(tmp > 0.1);
                angleOUT(slice,:) = [tOUT(itmp(1)+1),tOUT(itmp(1))];
            else
                angleOUT(slice,:) = [tOUT(1),tOUT(end)];
            end
        end
                
        %% Z
%         if selectedSlices(2)==selectedSlices(1)
%             tmp = ( DCMi.USliceLocation(slice) - DCMi.USliceLocation(selectedSlices(1)) );
%         else
        if selectedSlices(1) < 1
            tmpM = DCMi.USliceLocation(1) + abs(DCMi.spacing(3)) * (selectedSlices(1)-1);
        else
            tmpM = DCMi.USliceLocation(selectedSlices(1));
        end
        if selectedSlices(2) > length(DCMi.USliceLocation)
            tmpP = DCMi.USliceLocation(length(DCMi.USliceLocation)) + abs(DCMi.spacing(3)) * (selectedSlices(2)-length(DCMi.USliceLocation));
        else
            tmpP = DCMi.USliceLocation(selectedSlices(2));
        end
            
        if tmpP == tmpM
            tmp = 0;
        else
%             tmp = ( DCMi.USliceLocation(slice) - DCMi.USliceLocation(selectedSlices(1)) ) / ( DCMi.USliceLocation(selectedSlices(2)) - DCMi.USliceLocation(selectedSlices(1)) ); %%% may be higher than 1
            tmp = ( DCMi.USliceLocation(slice) - tmpM ) / ( tmpP - tmpM ); %%% may be higher than 1
        end
%         end
        I1 = NaN(Nr*samplingFactor,Nc*samplingFactor);
        I1(maskMYO==1) = tmp;
        I1(isnan(Inew)) = NaN;
        coords(:,:,slice,3) = I1;
        coordsList(slice) = tmp;
        
        if(slice == sliceToDraw)
            
            for d=1:3
                h2 = subplot(2,2,d+1); hold on;
                Itmp = squeeze(coords(:,:,slice,d));
                Ialpha = double(~isnan(Itmp));
                colormap(h2,jet(jetResolution));
                imagesc(Itmp,'XData',tmpX,'YData',tmpY,'AlphaData',Ialpha);
                axis equal; grid on;
                axis ij;
                zoom(1); %%% zoom(3) if not cropped
                axis manual;
                caxis([0,1]);
                if d==3
                    title(['Coord z = ',num2str(mean( Itmp(~isnan(Itmp)) ))]);
                end
            end

            drawnow;   
            
            fileNameFIG = [fileData.folderRoot,'P_DATA/0_PNG/',fileData.fileName(1:end-4),'_',num2str(slice),'.png'];
            saveas(gcf,fileNameFIG);
        end
    end
    close all;
end

coords = coords + 1; %%% between 1 and 2
coords(isnan(coords(:))) = 0;

data_IMG = data_IMG + 1; %%% between 1 and ...
data_IMG(isnan(data_IMG(:))) = 0;

data_SEG = data_SEG + 1; %%% between 1 and ...
data_SEG(isnan(data_SEG(:))) = 0;

save(fileData.fileNameMAT,'coords');
save(fileData.fileNameMAT,'coordsList','-append');
save(fileData.fileNameMAT,'angleOUT','-append');
save(fileData.fileNameMAT,'data_IMG','-append');
save(fileData.fileNameMAT,'data_SEG','-append');

end

