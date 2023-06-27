function [data_img,data_myo,data_inf,data_mvo] = cvi2Matrix(CONTOURS,DCM)

data_img = cell(2,1);
data_myo = cell(2,1);
data_inf = cell(2,1);
data_mvo = cell(2,1);

for mi=1:length(CONTOURS)

    %%% Separate contours
    if ~isempty(DCM{mi})
        slicesData = DCM{mi}.data;
        [Nr,Nc,Nz] = size(slicesData);

        data_img{mi} = zeros(Nr,Nc,Nz);
        data_myo{mi} = zeros(Nr,Nc,Nz);
        data_inf{mi} = zeros(Nr,Nc,Nz);
        data_mvo{mi} = zeros(Nr,Nc,Nz);

        for slice=1:Nz

            slicePos = slice;

            dataEpi    = CONTOURS{mi}.data{slicePos,1}; %%% xy
            dataEndo   = CONTOURS{mi}.data{slicePos,2};
            dataEnh    = CONTOURS{mi}.data{slicePos,4};
            dataReflow = CONTOURS{mi}.data{slicePos,5};

            data_img{mi}(:,:,slice) = squeeze(slicesData(:,:,slicePos));

            if slice <= size(CONTOURS{mi}.data,1)
                %%% Myocardium
                if(~isempty(dataEpi))
                    x1 = (dataEpi{1}(:,1) - DCM{mi}.offset(1))/DCM{mi}.spacing(1) + 1;
                    y1 = (dataEpi{1}(:,2) - DCM{mi}.offset(2))/DCM{mi}.spacing(2) + 1;
                    maskEpi = poly2mask(x1,y1,Nr,Nc);

                    if ~isempty(dataEndo)
                        x1 = (dataEndo{1}(:,1) - DCM{mi}.offset(1))/DCM{mi}.spacing(1) + 1;
                        y1 = (dataEndo{1}(:,2) - DCM{mi}.offset(2))/DCM{mi}.spacing(2) + 1;
                        maskEndo = poly2mask(x1,y1,Nr,Nc);
                        data_myo{mi}(:,:,slice) = maskEpi.*(~maskEndo);
                    else
                        data_myo{mi}(:,:,slice) = maskEpi;
                    end


                end
                %%% Infarct + MVO
                maskEnh = zeros(Nr,Nc);
                maskReflow = zeros(Nr,Nc);
                if (~isempty(dataEpi))
                    for sbi=1:length(dataEnh)
                        x1 = (dataEnh{sbi}(:,1) - DCM{mi}.offset(1))/DCM{mi}.spacing(1) + 1;
                        y1 = (dataEnh{sbi}(:,2) - DCM{mi}.offset(2))/DCM{mi}.spacing(2) + 1;
                        maskEnh = maskEnh + poly2mask(x1,y1,Nr,Nc);
                    end
                    maskEnh = maskEnh>0 ;

                    for sbi=1:length(dataReflow)
                        x1 = (dataReflow{sbi}(:,1) - DCM{mi}.offset(1))/DCM{mi}.spacing(1) + 1;
                        y1 = (dataReflow{sbi}(:,2) - DCM{mi}.offset(2))/DCM{mi}.spacing(2) + 1;
                        maskReflow = maskReflow + poly2mask(x1,y1,Nr,Nc);
                    end
                    maskReflow = maskReflow>0 ;

                    data_inf{mi}(:,:,slice) = maskEnh .* (~maskReflow) .* data_myo{mi}(:,:,slice);
                    data_mvo{mi}(:,:,slice) = maskReflow .* data_myo{mi}(:,:,slice);
                end
                                
            end
        end
    end
end

end
