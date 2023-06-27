function plotContours(CONTOURS,DCM)

typeList = {'MVO','LATE'};
colorList = {[0,1,0],[1,0,0],[0,0,1],[1,0,1],[0,1,1],[0.5,0.5,0.5]};
%%% 1=epi / 2=endo / 3=ROImyo / 4=ROIlesion / 5=lesion / 6=?

figure('units','normalized','position',[0,0,0.25,1],'color','w');
for mi=1:2 %%% MVO / LATE
    if(~isempty(CONTOURS{mi}))
        
        subplot(2,1,mi); hold on;
        Nz = size(DCM{mi}.data,3);
        for sci=1:Nz
            
            for ci=1:size(CONTOURS{mi}.data,2)
                if sum(ci==[4,5])>0
                    lineT = 3;
                else
                    lineT = 1;
                end
                for sbi=1:length(CONTOURS{mi}.data{sci,ci}) %%% if several contours for a given structure (ex: MVO)
                    tmpP = CONTOURS{mi}.data{sci,ci}{sbi};
                    if(~isempty(tmpP))
                        plot3(tmpP(:,1),tmpP(:,2),tmpP(:,3) ,'Color',colorList{ci},'LineWidth',lineT);
                    end
                end
            end
        end
        axis equal; % axis([MIN(1),MAX(1),MIN(2),MAX(2),MIN(3),MAX(3)]);
        view(20,5); grid on;
        title(typeList{mi});
    end
end

for mi=1:2 %%% MVO / LATE
    if(~isempty(CONTOURS{mi}))
        figure('units','normalized','position',[0.25,0.5*(mi-1),0.75,0.45],'Name',[typeList{mi},' - ',DCM{mi}.folder],'color','w');
        colormap(gray);
        Np = floor( sqrt(size(DCM{mi}.data,3)) );
        [Nr,Nc,Nz] = size(DCM{mi}.data);
        for sci=1:Nz
            
            slicePos = sci;
                        
            h = subplot(Np+1,Np+1,Nz-sci+1); hold on;

            tmpI = squeeze( DCM{mi}.data(:,:,slicePos) );
            tmpX = ([1,Nc]-1)*DCM{mi}.spacing(1) + DCM{mi}.offset(1);
            tmpY = ([1,Nr]-1)*DCM{mi}.spacing(2) + DCM{mi}.offset(2);

            imagesc(tmpI,'XData',tmpX,'YData',tmpY);
            
            %%% epi ROI for brightness display
            dataEpi = CONTOURS{mi}.data{slicePos,1};
            if(~isempty(dataEpi))            
                x1 = (dataEpi{1}(:,1) - DCM{mi}.offset(1))/DCM{mi}.spacing(1) + 1;
                y1 = (dataEpi{1}(:,2) - DCM{mi}.offset(2))/DCM{mi}.spacing(2) + 1;
                maskEpi = poly2mask(x1,y1,Nr,Nc);
%                 set(h,'CLim',[min(tmpI(maskEpi))*0.5,max(tmpI(maskEpi))*1.5]);
                set(h,'CLim',[min(tmpI(maskEpi))*0.9,max(tmpI(maskEpi))*1.1]);
            end
            
            for ci=1:size(CONTOURS{mi}.data,2)
                for sbi=1:length(CONTOURS{mi}.data{slicePos,ci}) %%% if several contours for a given structure (ex: MVO)
                    tmpP = CONTOURS{mi}.data{slicePos,ci}{sbi}-2; %-1;
                    if(~isempty(tmpP))
                        plot(tmpP(:,1),tmpP(:,2),'Color',colorList{ci});
                    end
                end
            end
            axis equal; grid on; axis ij;
            axis off;
            zoom(1); %%% zoom(3) if not cropped
            title(['Slice ',num2str(slicePos),' / ',num2str(DCM{mi}.USliceLocation(sci))],'FontSize',6);
        end
    end
end

end
