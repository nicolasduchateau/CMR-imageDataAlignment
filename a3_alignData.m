clear all;
close all;
clc;

addpath('./Functions');

folderRoot = '../DATA/';

typeList = {'EGE','LGE'};

singleSlice = 0;

%% Create reference anatomy

if singleSlice==0
    fileNameREF = [folderRoot,'P_DATA/COORDS/','CoordsAndData_REF.mat'];
else
    fileNameREF = [folderRoot,'P_DATA/COORDS/','SLICE_CoordsAndData_REF.mat'];
end

if exist(fileNameREF,'file') ~= 2
    Nz = 21;
    zList = linspace(0,0.9,Nz);
    zList = fliplr(zList);

    if singleSlice==1
        zList = repmat(zList( (Nz-1)/2+1 ),[1,Nz]); %%% SINGLE SLICE
    end

    Nr = 80;
    Nc = 80;
    Nz = length(zList);
    DCMi.data = NaN(Nr,Nc,Nz);
    DCMi.offset = [0,0,0];
    DCMi.spacing = [1,1,5]';
    DCMi.USliceLocation = ((1:Nz)-1)*DCMi.spacing(3);

    Np = 200;
    CONTOURSi.data = cell(Nz,6);
    t = linspace(0,2*pi,Np+1)' + 4*pi/3; t(end) = [];
    b = 1;
    a = 0.5;
    rList = ( a * (1-(zList/b).^2) ).^(1/2);
    originListi = zeros(Nz,2);
    for zi=1:Nz
        z = zList(zi);
        radEPI = 50*rList(zi);
        radENDO = 30*rList(zi);

        CONTOURSi.data{zi,1}{1} = [radEPI*cos(t)+Nc/2,radEPI*sin(t)+Nc/2,repmat(z,[Np,1])];
        CONTOURSi.data{zi,2}{1} = [radENDO*cos(t)+Nr/2,radENDO*sin(t)+Nr/2,repmat(z,[Np,1])];
        originListi(zi,:) = CONTOURSi.data{zi,1}{1}(1,[2,1]);
    end
    
    fileData.folderRoot = folderRoot;
    if singleSlice==0
        fileData.fileName = ['CoordsAndData_REF.mat'];
    else
        fileData.fileName = ['SLICE_CoordsAndData_REF.mat'];
    end
    fileData.fileNameMAT = [fileData.folderRoot,'P_DATA/COORDS/',fileData.fileName];
    selectedSlices = [1,Nz];
    myoOpening = NaN(Nz,2,2);

    samplingFactor = 4;
    computeLocalCoords(originListi,myoOpening,selectedSlices,fileData,CONTOURSi,DCMi,samplingFactor);

    if singleSlice==0
        fileData.fileName = ['CoordsAndData_REF_LowRes.mat'];
    else
        fileData.fileName = ['SLICE_CoordsAndData_REF_LowRes.mat'];
    end
    fileData.fileNameMAT = [fileData.folderRoot,'P_DATA/COORDS/',fileData.fileName];
    samplingFactor = 1;
    computeLocalCoords(originListi,myoOpening,selectedSlices,fileData,CONTOURSi,DCMi,samplingFactor);
end

%% Resample each subject's data to the reference

Nfeat = 5; %%% 1=imageData / 2=binaryInfarct / 3-4-5=coords

load(fileNameREF);

coords = coords - 1;
coords(coords==-1) = NaN;

idxREF = find(~isnan(coords(:,:,:,3)));
tmp_coordsREF = NaN(length(idxREF),3);
for d=1:3
    tmp = squeeze(coords(:,:,:,d));
    tmp_coordsREF(:,d) = tmp(idxREF);
end            
[NrREF,NcREF,NzREF,~] = size(coords);                

for ki=1%1:140
    k = ki;
    disp(k);
    for mi=1:2 %%% 1=EGE / 2=LGE
        if singleSlice==0
            fileName = [folderRoot,'P_DATA/COORDS/','CoordsAndData_',typeList{mi},'_',sprintf('%03.0f',k),'.mat'];
            fileNameMAT = [folderRoot,'P_DATA/RESAMP/RESAMP_',typeList{mi},'_',sprintf('%03.0f',k),'.mat'];
        else
            fileName = [folderRoot,'P_DATA/COORDS/','CoordsAndData_',typeList{mi},'_',sprintf('%03.0f',k),'.mat'];
            fileNameMAT = [folderRoot,'P_DATA/RESAMP/SLICE_RESAMP_',typeList{mi},'_',sprintf('%03.0f',k),'.mat'];
        end

        if exist(fileName,'file')==2 && exist(fileNameMAT,'file')~=2
            load(fileName);

            coords = coords - 1;
            coords(coords==-1) = NaN;
            data_IMG = data_IMG - 1;
            data_IMG(data_IMG==-1) = NaN;
            data_SEG = data_SEG - 1;
            data_SEG(data_SEG==-1) = NaN;
            if mi==1
                data_SEG(data_SEG<2.5) = 1;
            end

            [Nr,Nc,Nz,~] = size(coords);                

            IN_data = cell(Nfeat,1);
            idx = find(~isnan(coords(:,:,:,1)));
            tmp_coords = NaN(length(idx),3);
            tmp = data_IMG(:,:,:);
            IN_data{1} = tmp(idx);
            tmp = data_SEG(:,:,:);
            IN_data{2} = tmp(idx);
            for d=1:3
                tmp = squeeze(coords(:,:,:,d));
                tmp_coords(:,d) = tmp(idx);
                IN_data{2+d} = tmp_coords(:,d);
            end

            INcoords = unique(tmp_coords(:,3));
            INcoordsREF = unique(tmp_coordsREF(:,3));

            OUT_data = cell(Nfeat,1);
            for feat=1:Nfeat
                OUT_data{feat} = NaN(size(tmp_coordsREF,1),1);
            end
            for feat=[1,2] %1:Nfeat
                F = scatteredInterpolant(tmp_coords(:,1),tmp_coords(:,2),tmp_coords(:,3),IN_data{feat}(:),'linear','nearest');
                vq = F(tmp_coordsREF(:,1),tmp_coordsREF(:,2),tmp_coordsREF(:,3));
                OUT_data{feat}(:) = vq;
            end
            for feat=[3,4,5] %1:Nfeat
                OUT_data{feat}(:) = tmp_coordsREF(:,feat-2);
            end

            %%%%% remove openings data
            ids = find(~isnan(angleOUT(:,1)));
            if(~isempty(ids))
                idsREF = find(INcoordsREF > coordsList(ids(1)-1)); %%% ids(1)-1 = slice before first slice with opening
                %%% interpolate angleOUT at these slices
                tmp = angleOUT(ids,:);
                idx2 = find(tmp(:,1) > 0.5);
                tmp(idx2,1) = tmp(idx2,1) - 1;
                if length(ids)==1
                    vq = repmat(tmp,[length(idsREF),1]);
                else
                    vq = NaN(length(idsREF),2);
                    for p=1:2
                        vq(:,p) = interp1(coordsList(ids)',tmp(:,p)',INcoordsREF(idsREF)','linear','extrap')';
                    end
                end
            end

            tmp = cell(2,1);
            tmp{1} = IN_data;
            tmp{2} = OUT_data;
            idT = {idx,idxREF};
            sizeT = {[Nr,Nc,Nz],[NrREF,NcREF,NzREF]};
            plotData = cell(2,1);
            for c=1:2
                plotData{c} = cell(Nfeat,1);
                for feat=1:Nfeat
                    plotData{c}{feat} = NaN(sizeT{c});
                    plotData{c}{feat}(idT{c}) = tmp{c}{feat};
                end
            end

            if(~isempty(ids))
                for c=2
                    tmpOUT = ones([NrREF,NcREF,NzREF]);
                    for si=1:length(idsREF)
                        slice = idsREF(si);
                        tmp = plotData{c}{4}(:,:,slice); %%% if CIRCcoord between the two anglesOUT = NaN 
                        if( vq(si,1) < 0 )
                            tmp = 1 - ( (tmp > (1+vq(si,1))) + (tmp < vq(si,2)) );
                        else
                            tmp = 1 - (tmp > vq(si,1)) .* (tmp < vq(si,2));
                        end
                        tmp(tmp==0) = NaN;
                        tmpOUT(:,:,slice) = tmp;
                    end
                    for feat=1:Nfeat
                        plotData{c}{feat} = plotData{c}{feat} .* tmpOUT;

                        %%% remove isolated points
                        tmpI = ~isnan(plotData{c}{feat});
                    end
                end
            end

            %%% Remove data outside extreme slides
            tmpOUT = ones([NrREF,NcREF,NzREF]);
            deltaS = INcoords(2)-INcoords(1);
            for si=1:NzREF
                if( INcoordsREF(si)<INcoords(1)-deltaS || INcoordsREF(si)>INcoords(end)+deltaS )
                    tmpOUT(:,:,si) = NaN;
                end
            end
            for c=2
                for feat=1:Nfeat
                    plotData{c}{feat} = plotData{c}{feat} .* tmpOUT;
                end
            end

            %%% Plotting
            for variableToPlot = 1:2
                figure('units','normalized','position',[0,0,1,1],'Visible','on');
                for c=1:2
                    if c==1
                        Cplot = coordsList;%INcoords;
                    else
                        Cplot = INcoordsREF;
                    end
                    for d=variableToPlot
                        for s=1:length(Cplot)
                            h1 = subplot(4,11,s + 22*(c-1)); hold on;
                                                                                    
                            Itmp = squeeze( plotData{c}{d}(:,:,s) );
                            Ialpha = double(~isnan(Itmp));
                            if d==1
                                colormap(h1,gray);
                                if( sum(~isnan(Itmp(:))) > 0 )
    %                                 caxis([min(Itmp(:))*0.5,max(Itmp(:))*1.5]);
                                    caxis([min(Itmp(:))*0,max(Itmp(:))*1]);
                                end
                            elseif d==2
                                colormap(h1,jet);
                                caxis([0,3]);
                            else
                                colormap(h1,jet);
                                caxis([0,1]);
                            end
                            imagesc(Itmp,'AlphaData',Ialpha);
                            B = bwboundaries(Ialpha);

                            if ~isempty(B)
                                plot(B{1}(:,2),B{1}(:,1),'r');
                                if length(B) == 2
                                    plot(B{2}(:,2),B{2}(:,1),'r');
                                end
                            end

                            
                            tmp3 = plotData{c}{3}(:,:,s);
                            tmp4 = plotData{c}{4}(:,:,s);
                            idx = find( (tmp3(:)>0.8) .* (tmp4(:)<0.01) );
                            [Ir,Jr] = ind2sub([size(plotData{c}{3},1),size(plotData{c}{3},2)],idx);
                            Ir = mean(Ir);
                            Jr = mean(Jr);
                            plot(Jr,Ir,'.g','MarkerSize',20);

                            title(num2str(Cplot(s)));
                            axis equal; grid on; axis ij;
                            set(gca,'FontSize',6);
                            zoom(1);
                        end
                    end
                end
                fileNameFIG = [folderRoot,'P_DATA/1_PNG/R_',num2str(k),'_',typeList{mi},'_',num2str(variableToPlot),'.png'];
                saveas(gcf,fileNameFIG);                
                close all;

            end

            ResampledData = plotData{2};
            save(fileNameMAT,'ResampledData');            
        end
    end
end




