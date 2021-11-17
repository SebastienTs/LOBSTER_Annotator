clear all;
close all;
clc;
addpath(genpath('functions'));
global ImgSize;
global RGB;
global LblDil;

ImgsPath = [uigetdir('Images/','Select 2D TIFF images folder') '\'];

%% Configuration
ModelName = 'Model';
ModelsPath = [ImgsPath '\Models'];
Train       = 1;
Pred        = 0;

%% Model parameters
ImgSize = 512;
RGB = 0;
LblDil = 0;
classNames = ["a","b"];
labelIDs = [0 255];
NumFirstLayerChan = 16;
EncDepth = 6;
Loss = 'gdice'; %% 'bce', 'wbce' or 'gdice'
w_b = 0.5;
w_o = 0.5;

%% Training parameters
Epochs = 200;
NVal = 2;
MiniBatch = 4;

%% Configure data store
imdsAll = imageDatastore(fullfile(ImgsPath,'*.tif'));
NImgs = numel(imdsAll.Files);
pxds = pixelLabelDatastore(fullfile(ImgsPath,'*_ref.png'),classNames,labelIDs);
wrapper = @(x) strrep(x,'_ref.png','.tif');
imds = imageDatastore(fullfile(ImgsPath,'*.tif'));
imds.Files = cellfun(wrapper,pxds.Files,'UniformOutput',false);
NRef = numel(imds.Files);

%% Dialog box
defaultValues = {num2str(ImgSize),num2str(Epochs),num2str(NVal),num2str(MiniBatch),'1'};
titleBar = 'UNET_Train';
userPrompt = {'Image size (128-1024)', 'Epochs (5-1000)',['Number of validation images 1-' num2str(NImgs-1)],'MiniBatch (1-16)'};
caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
if ~isempty(caUserInput)
    ImgSize = boundvar(str2double(caUserInput{1}),128,1024,512);
    Epochs = boundvar(str2double(caUserInput{2}),5,1000,200);
    NVal = boundvar(str2double(caUserInput{3}),1,NImgs-1,1);
    MiniBatch = boundvar(str2double(caUserInput{4}),1,64,1);
end
ValFreq = 10*round((NRef-NVal)/MiniBatch);

%% Training
if Train
    
    %% Label + Image DataStore
    augmenter = imageDataAugmenter('RandXReflection',true,'RandYReflection',true);
    ds = pixelLabelImageDatastore(imds,pxds,'DataAugmentation',augmenter);
    dstrain = partitionByIndex(ds,1:NRef-NVal);
    dsval = partitionByIndex(ds,NRef);
    dstrain = transform(dstrain,@commonPreprocessing);
    dsval = transform(dsval,@commonPreprocessing);
    
    %% Create U-NET
    if RGB == 0
        lgraph = unetLayers([ImgSize ImgSize],numel(labelIDs),'FilterSize',[3 3],'EncoderDepth',EncDepth,'NumFirstEncoderFilters',NumFirstLayerChan);
    else
        lgraph = unetLayers([ImgSize ImgSize 3],numel(labelIDs),'FilterSize',[3 3],'EncoderDepth',EncDepth,'NumFirstEncoderFilters',NumFirstLayerChan);
    end
    
    %% Set network loss function for training
    switch Loss
        case 'gdice'
            %% Generalized DICE (class weighting)
            layer = dicePixelClassificationLayer('Name','Dice-Layer');
            lgraph = replaceLayer(lgraph,'Segmentation-Layer',layer);
        case 'wbce'
            %% Binary cross-entropy + class weighting
            tbl = countEachLabel(pxds);
            layer = pixelClassificationLayer('Name','Segmentation-Layer','Classes',tbl.Name,'ClassWeights',[w_b w_o]);
            lgraph = replaceLayer(lgraph,'Segmentation-Layer',layer);
        otherwise
           %% Binary cross-entropy 
    end     
     
    %% Train U-NET
    options = trainingOptions('adam','InitialLearnRate',5e-4,'ValidationData',dsval,'MiniBatchSize',MiniBatch,'MaxEpochs',Epochs,'ValidationFrequency',ValFreq,'VerboseFrequency',ValFreq,'Shuffle','never');
    net = trainNetwork(dstrain,lgraph,options);

    %% Save net
    save([ModelsPath '\' ModelName],'net');    
   
end

%% Prediction & evaluation
if Pred
    
    %% Load net
    load([ModelsPath '\' ModelName]);

    %% Predict
    dstest = transform(imdsAll,@commonPreprocessing);
    Pred = predict(net,dstest);
    P = Pred(:,:,2,:);
    
    %% Save images
    for i=1:size(P,4)
        imwrite(uint8(255*(P(:,:,i)>0.5)),strrep(imdsAll.Files{i},'.tif','_pre.png'));
    end
    
end

function dataOut = commonPreprocessing(data)
    global ImgSize;
    global RGB;
    global LblDil;
    dataOut = cell(size(data));
    if size(data,2) == 2
        for idx = 1:size(data,1)
            %% Images
            temp = data{idx,1};
            temp = temp{1};
            if size(temp,3)>1 && RGB==0
                temp = mean(temp,3);
            end
            temp = imresize(temp,[ImgSize,ImgSize],'bilinear');
            dataOut{idx,1} = temp;
            %% Label masks
            temp = data{idx,2};
            temp = temp{1};
            if abs(LblDil)>1
                C = categories(temp);
                temp = double(temp);
                if LblDil>1
                    temp = imdilate(temp,true(LblDil));
                elseif LblDil<-1
                    temp = imerode(temp,true(abs(LblDil)));
                end
                temp = categorical(temp,[1 2],{C{1},C{2}});
            end
            temp = imresize(temp,[ImgSize,ImgSize],'nearest');
            dataOut{idx,2} = temp;
        end
    else
        %% Images
        temp = data;
        if size(temp,3)>1 && RGB==0
            temp = mean(temp,3);
        end
        temp = imresize(temp,[ImgSize,ImgSize],'bilinear');
        dataOut = temp;
    end
end