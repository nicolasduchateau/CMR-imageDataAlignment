function [origine,myoOpening] = initializeOrigin(data_img,data_myo,data_inf,data_mvo,fileData,CONTOURSi,DCMi,selectedSlices)

[Nr,Nc,Nz] = size(data_img);
origine = zeros(Nz,2);
myoOpening = NaN(Nz,2,2);

displayAllSlices(data_img,data_myo,data_inf,data_mvo,CONTOURSi,DCMi);

h0 = figure('units','normalized','position',[0 0 0.2 0.4]); hold on;
slicesSubset = selectedSlices(2):-1:selectedSlices(1);
slicesSubset(slicesSubset<1) = [];
slicesSubset(slicesSubset>Nz) = [];
for slice=slicesSubset
    
    hold on; colormap gray;
    
    imagesc(data_img(:,:,slice));

    %%% epi ROI for brightness display
    dataEpi = CONTOURSi.data{slice,1}; %%% xy
    if(~isempty(dataEpi)) 
        tmpI = squeeze( data_img(:,:,slice) );
        x1 = (dataEpi{1}(:,1) - DCMi.offset(1))/DCMi.spacing(1) + 1;
        y1 = (dataEpi{1}(:,2) - DCMi.offset(2))/DCMi.spacing(2) + 1;
        maskEpi = poly2mask(x1,y1,Nr,Nc);
%         caxis([min(tmpI(maskEpi))*0.5,max(tmpI(maskEpi))*1.5]);
        caxis([min(tmpI(maskEpi))*0.9,max(tmpI(maskEpi))*1.1]);
    end
    
    BW =  data_myo(:,:,slice);
    [B,~] = bwboundaries(BW,'holes');
    
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2);
        axis equal;
        axis ij;
    end
    
    x = 0;
    y = 0;
    if ~isempty(B)
        isOK = 0;
        title('Draw 1 point RV-LV (+ 2 points myo. opening)');
        while isOK==0
            [x,y] = getpts; %%% xy but in Nc,Nr coords = indices
            if sum( length(x)==[1,3] ) > 0
                isOK = 1;
            end
        end
    end
    origine(slice,:) = [y(1),x(1)]; %%% BEWARE, indices stored
    if length(y) == 3
        myoOpening(slice,:,:) = [y(2:3),x(2:3)]; %%% BEWARE, indices stored
    end
    
end

ind =  find(origine(:,1) < 0); %%% force not defining origin on some slices
origine(ind,:) = repmat([0,0],[length(ind),1]);

close(h0);

if( length(slicesSubset) > 1)
    maxSlice = max(find(origine(:,1)> 0));
    minSlice = min(find(origine(:,1)> 0));
    middleSlice = minSlice + round((maxSlice - minSlice +1) / 2 );

    h1 = figure('units','normalized','position',[0 0 0.6 1]); hold on;

    imagesc(data_img(:,:,middleSlice));
    colormap gray;
    axis equal; axis ij;

    title('Draw RV to LV cut (1 line)');
    isOK = 0;
    while isOK==0
        h = imline;
        pos = h.getPosition;
        if( (size(pos,1)==2) && (pos(1,1)~=pos(2,1)) ) %%% numPoints
            isOK = 1;
        end
    end

    vectorCoef = [pos(2,1) - pos(1,1) ,pos(2,2) - pos(1,2)]; %%% xy in pixels

    distances = 0:-0.1:-10;
    i = vectorCoef(2) * distances + pos(1,2);
    j = vectorCoef(1) * distances + pos(1,1);
    ind = find(i >=Nr | i < 0 | j >=Nc | j < 0,1,'first');
    p1 = [i(ind-1),j(ind-1)];
    distances = 0:0.1:10;
    i = vectorCoef(2) * distances + pos(1,2);
    j = vectorCoef(1) * distances + pos(1,1);
    ind = find(i >=Nr | i < 0 | j >=Nc | j < 0,1,'first');
    p2 = [i(ind-1),j(ind-1)]; %%% ij

    plot(p1(2),p1(1),'r*'); %%% xy
    plot(p2(2),p2(1),'r*');

    [indDrawline,~] = drawline(p1,p2,[Nr,Nc]);

    chamberView = zeros(Nz,length(indDrawline));
    for slice=Nz:-1:1
        data = data_img(:,:,slice);
        chamberView(slice,:) = data(indDrawline);    
    end

    close(h1);

    h2 = figure('units','normalized','position',[0 0 0.3 0.3]); hold on;
    imagesc(chamberView);
    axis([1,size(chamberView,2),1,size(chamberView,1)]);
    pbaspect([size(chamberView,2)*DCMi.spacing(1) , size(chamberView,1)*abs(DCMi.spacing(3)) , 1]); %%% spacing aspect ratio
    axis off;
    % zoom(2);

    set(gca,'YDir','normal')
    colormap gray;
    SegEndoEpi = find(origine(:,1)> 0);
    hold on;
    plot((max(SegEndoEpi))*ones(length(indDrawline),1),'r', 'LineWidth', 2);
    plot((min(SegEndoEpi))*ones(length(indDrawline),1),'r', 'LineWidth', 2);

    plot((selectedSlices(2))*ones(length(indDrawline),1),'y--', 'LineWidth', 2);
    plot((selectedSlices(1))*ones(length(indDrawline),1),'y--', 'LineWidth', 2);
    title('PAUSED - press key to continue');
    pause;

    if selectedSlices(1) < minSlice
        tmp = max(1,selectedSlices(1));
        origine(tmp:minSlice-1,:) = repmat(origine(minSlice,:),[minSlice-tmp,1]);
    else
        origine(minSlice:selectedSlices(1)-1,:) = repmat([0,0],[selectedSlices(1)-minSlice,1]);
    end
    if selectedSlices(2) > maxSlice
        tmp = min(Nz,selectedSlices(2));
        origine(maxSlice+1:tmp,:) = repmat(origine(maxSlice,:),[tmp-maxSlice,1]);
    else
        origine(selectedSlices(2)+1:maxSlice,:) = repmat([0,0],[maxSlice-selectedSlices(2),1]);
    end
end

save(fileData.fileName,'origine');
save(fileData.fileName,'myoOpening','-append');

close all;

end
