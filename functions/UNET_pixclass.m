function UNET_pixclass(ImgsPath,ImgSz,ClassBck,ClassObj)

    global ImgSize;
    global RGB;

    %% Configuration
    Files = dir(strcat([ImgsPath '/Models/*.mat']));
    ModelName = Files(1).name;
    
    %% Model parameters
    ImgSize = [512 512];
    RGB = 0;

    %% Configure data store
    imdsAll = imageDatastore(fullfile(ImgsPath,'*.tif'));

    %% Load net
    load([ImgsPath '\Models\' ModelName],'net');
    disp(['Classifying pixels with model ' ModelName]);
    
    %% Predict
    dstest = transform(imdsAll,@commonPreprocessing);
    Pred = predict(net,dstest);
    P = Pred(:,:,2,:);

    %% Save images
    for i=1:size(P,4)
        Img = uint8(double(ClassObj)*(P(:,:,i)>=0.5)+double(ClassBck)*(P(:,:,i)<0.5));
        Img = imresize(Img,ImgSz,'nearest');
        imwrite(Img,strrep(imdsAll.Files{i},'.tif','.png'));
    end
    
end

function dataOut = commonPreprocessing(data)
    global ImgSize;
    global RGB;
    temp = data;
    if size(temp,3)>1 && RGB==0
        temp = mean(temp,3);
    end
    temp = imresize(temp,ImgSize,'bilinear');
    dataOut = temp;
end