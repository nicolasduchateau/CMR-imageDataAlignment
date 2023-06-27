function displayAllSlices(data_img,data_myo,data_inf,data_mvo,CONTOURSi,DCMi)
%% Display on 1 figure all the slices from a case with myocardium and infarct zone
% and in an other the histogram of the slices

[Nr,Nc,Nz] = size(data_img);

hd = figure('units','normalized','position',[0.6 0 0.4 1]);
% if ~isempty(originList)
%     ori = originList.origine;
% end

% if ~exist('selectedSlice','var')
    % third parameter does not exist, so default it to something
    selectedSlice = [1,Nz];
% end

n = ceil((selectedSlice(2)-selectedSlice(1)+1)/5);

width = 1/5-0.01;
height = 1/n-0.01;

for i=selectedSlice(1):1:selectedSlice(2)
    
    colormap gray;
    posI = i- selectedSlice(1) +1;
    col = width * (mod(posI-1,5));
    row = 1 - height * ceil(posI/5);
    pos1 = [col row width height];
    
    h = subplot('Position',pos1); hold on;
    
    title(num2str(posI),'Position',[pos1(1),pos1(2)]);
   
    imagesc(data_img(:,:,i));
    axis equal; axis off; axis ij;
    
    %%% Myocardium
    BW = imbinarize(data_myo(:,:,i));
    [B,~] = bwboundaries(BW,'holes');
    
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1);
    end
    
    %%% epi ROI for brightness display
    dataEpi = CONTOURSi.data{i,1};
    if(~isempty(dataEpi)) 
        tmpI = squeeze( data_img(:,:,i) );
        x1 = (dataEpi{1}(:,1) - DCMi.offset(1))/DCMi.spacing(1) + 1;
        y1 = (dataEpi{1}(:,2) - DCMi.offset(2))/DCMi.spacing(2) + 1;
        maskEpi = poly2mask(x1,y1,Nr,Nc);
        set(h,'CLim',[min(tmpI(maskEpi))*0.5,max(tmpI(maskEpi))*1.5]);
    end
    
    idx = [];
    %%% Infarct
    if ~isempty(data_inf)
        idx = find((data_inf(:,:,i) == 1));
    end
    inf = size(idx);
    [I,J] = ind2sub([Nr,Nc],idx);
    s = size(I);
    if(s(1)> 1 ) %%% Check if the myocardium is segmented or not
        plot(J,I,'.r');
    end
    indInf = idx;
    
    idx = [];
    %%% MVO
    if ~isempty(data_inf)
        idx = find((data_mvo(:,:,i) == 1));
    end
    inf = size(idx);
    [I,J] = ind2sub([Nr,Nc],idx);
    s = size(I);
    if(s(1)> 1 ) %%% Check if the myocardium is segmented or not
        plot(J,I,'.g');
    end
    
%     %%% Plot origin
%     if ~isempty(originList)
% %         plot(ori(i,1),ori(i,2),'g+','Markersize',5);
%     end
    indMyo = find((data_myo(:,:,i)==1));
        
    zoom(1);
end

drawnow;

end
