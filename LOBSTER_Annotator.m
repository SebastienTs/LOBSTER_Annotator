clear all;
close all;
clc;
%addpath(genpath('functions'));

%% ML classifier configuration
BlckSizeXY = 3;         %% XY block size (3D only)
BlckSizeZ = 3;          %% Z block size (3D only)
Scales = [1 3 6];       %% Scale, or scale geometrical progression (3D only)
NShifts = 1;            %% Number of directional shifts at each scale
predchunks = 250000;    %% Prediction chunk size (number of blocks)
NTrees = 25;            %% Random forest number of trees (only for >2 classes)
BckThr = 0;             %% Do not classify blocks with mean intensity below this threhold
ThrLvl = 0.5;           %% Voxels below block local threshold are set to background (use negative values for bright background)
SpaFeat = 0;            %% Use spatial features

%% Process single image, folder of 2D images or 3D image folder?
choice = menu('LOBSTER Annotator','3D TIFF file (4 GB max)','TIFF series (3D image folder)','2D TIFF images (folder)');
switch choice
    case 1
        FolderMode = 0;
        ImageSeries = 0;
    case 2
        FolderMode = 0;
        ImageSeries = 1;
    case 3
        FolderMode = 1;
        ImageSeries = 0; 
end

%% FolderMode configuration
if FolderMode == 1
    selpath = [uigetdir('Images/','Select 2D TIFF images folder') '\'];
    Files = dir(strcat(selpath,'*.tif'));
    NImages = numel(Files);
    disp(['Found ' num2str(NImages) ' images']);
    wanddim = 0;            %% Default to 2D wand
    BlckSizeXY = 2;         %% No padding
    BlckSizeZ = 2;          %% No padding
    Scales = [1 2 4 8 16]; %% Image assumed multiple of 16  
else      
    wanddim = 1;        %% Default to 3D wand
    NImages = 1;        %% Only 1 image to process
end

%% Settings
LUT = uint8([0 50 75 100 125 200 220 250]); %% Classes intensity levels
global LastColor;
LastColor = LUT(2);     %% Current class for brushes
Alpha = 1;              %% Annotation mask opacity
SizeLimit = 2^26;       %% Recommended maximum target image size (64 MB)
PitchLimit = 2^20;      %% Recommended maximum XY size (1024 x 1024)
NCols = 16;             %% Number of colors for 3D rendering LUT
PointSize = 5;          %% 3D rendering point size
ZRatio = 2;             %% Default Z Ratio
projmode = 0;           %% Default: No local Z projection
RunProj = 15;           %% Default local Z projection depth
ZSpan = 256;            %% Default: render all Z slices
Inv = 0;                %% Video inversion

%% Default values (can be changed from parameters dialog boxes)
Sig = 4;                %% Blobs: default Gaussian blur radius
LocalMaxBox = 7;        %% Blobs: defualt local maxima search box
MinLoG = 0.1;           %% Blobs: default sensitivity
FillMode = 0;           %% Blobs: default growing algorithm (0:watershed 1:propagate)
wandblur = 1;           %% Wand: default Gaussian blur radius
wandtol = 0.5;          %% Wand: default tolerance
wanddst = Inf;          %% Wand: default propagation distance
%Tracedepth = 10;
wandmode = 1;           %% Wand: default thresholding mode (0:low/up, 1:low)
splitmode = 1;          %% Wand: default splitting algorithm (0: regiondescent, 1: 3D watershed)
gaussraddstmap = 2;     %% Wand: default Gaussian radius of distance map for splitting
splitdst = 100;         %% Wand: default regiondescent splitting distance
growdst = 100;          %% Wand: default regiondescent growing distance
PreCloseRad = 1;        %% Skeleton: default morphological closing radius
Min2DHolesArea = 15;    %% Skeleton: default minimum 2D area closing 
MinVol = 50;            %% Skeleton: default minimum isolated skeleton volume
MinBrchLgth = 7;        %% Skeleton: default minimum branch length
PruneIter = 1;          %% Skeleton: default number of pruning iterations
Ignore4Way = 0;         %% Skeleton: default ignore 4-way branches
global RunProjViewer;
RunProjViewer = 0;      %% Current depth of running projection
global Polarity;
Polarity = 1;           %% Smart brush bright/dark objects
FastSklLbl = 0;         %% Default skeleton labeling
global KeyRead;         %% Used to read currently pressed key
AnnSource = 0;          %% Use regular or reference annotation file names

%% Initialization
StepXY = floor(BlckSizeXY/2);
StepZ = floor(BlckSizeZ/2);
ShowAnnot = true;
MaskMode = false;
brchmode = 1;
predinterp = 1;

%% Compute LUT3d (3D rendering)
LUT3d = round(linspecer(NCols)*255);
rng(0);
LUT3d = LUT3d(randperm(NCols),:);

%% 6 connectivity structure element
se = uint8(zeros(3,3,3));
se(:,:,1) = [0 0 0;0 1 0;0 0 0];
se(:,:,2) = [0 1 0;1 1 1;0 1 0];
se(:,:,3) = [0 0 0;0 1 0;0 0 0];

%% Image filename filter
if ImageSeries == 1                     
    defaultValues = {''};
    titleBar = 'Filename filter';
    userPrompt = {'Optional filename substring: '};
    caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
    if ~isempty(caUserInput)
        Filter = caUserInput{1};
    end
end

%% Initialize figure window
handle = figure('units','normalized','outerposition',[0 0 1 1]);
set(handle, 'KeyPressFcn',@keyPress);
set(handle, 'MenuBar', 'none');
set(handle, 'ToolBar', 'none');
set(handle,'Color','k');

for it = 1:NImages
  
    if FolderMode == 0
        %% Image, existing annotations & display
        if ImageSeries == 0
            [file,lastpath] = uigetfile('Images/*.tif','Select 2D/3D TIFF file');
            fname = [lastpath file];
            info = imfinfo(fname);
            num_slices = numel(info);
            if num_slices == 1
                aname = [lastpath file(1:end-4) '_ann.png'];
            else
                aname = [lastpath file(1:end-4) '_ann.tif']; 
            end
        else
            flname = uigetdir('Images/','Select a 3D TIFF image folder');
            [filepath,name,ext] = fileparts(flname);
            aname = [filepath '\' name '_ann.tif'];
        end
    else
        fname = [selpath Files(it).name];
    end
    
    %% Display image index
    disp(['Image ' num2str(it) ' of ' num2str(NImages)]);    
    set(handle, 'Name',['Image ' num2str(it) ' of ' num2str(NImages)]);
    
    %% Load image
    disp('Loading image...');
    h = waitbar(0,'Loading image');
    if ImageSeries == 0
        info = imfinfo(fname);
        num_slices = numel(info);
        if num_slices == 1
            Mode2D = 1;
            BlckSizeZ = 1;
        else
            Mode2D = 0;
        end
        if FolderMode == 1 && num_slices > 1
            error('Folder mode is incompatible with 3D images');
        end
        OrigWidth = info(1).Width;
        OrigHeight = info(1).Height;
        
        %% Downsampling
        if FolderMode == 0
            if OrigWidth*OrigHeight*num_slices > SizeLimit 
                RecDownSamp = ceil(sqrt(OrigWidth*OrigHeight*num_slices/SizeLimit));
                str = inputdlg(sprintf('Processing might be slow for images exceeding 64 MB\nRecommended XY downsampling factor'),'Large image detected',1,{num2str(RecDownSamp)});
                DownSamp = max([str2double(str{1}) 1]);
            elseif OrigWidth*OrigHeight > PitchLimit
                RecDownSamp = ceil(sqrt(OrigWidth*OrigHeight/PitchLimit));
                str = inputdlg(sprintf('Interface might be slow for image slices exceeding 1 MB\nRecommended XY downsampling factor'),'Large image detected',1,{num2str(RecDownSamp)});
                DownSamp = max([str2double(str{1}) 1]);
            else
                DownSamp = 1;
            end
        else
            DownSamp = 1;
        end
        
        %% Initialize arrays
        Width = floor(info(1).Width/DownSamp);
        Height = floor(info(1).Height/DownSamp);
        if info(1).BitDepth == 8
            I = uint8(zeros(Height,Width,num_slices));
        else
            I = uint16(zeros(Height,Width,num_slices));
        end
        
        %% Loading images
        for kf = 1:num_slices
            h = waitbar(kf/num_slices);
            slice = imread(fname, kf, 'PixelRegion',{[1+floor(DownSamp/2),DownSamp,Height*DownSamp],[1+floor(DownSamp/2),DownSamp,Width*DownSamp]});
            %% Color mode
            if size(slice,3) > 1
                disp('Color images not supported, loading RGB average');
                slice = round(mean(slice,3));
            end
            I(:,:,kf) = slice;
        end
   
    else %% 3D image as series
        
        select = 1;
        while select
            %% Analyze image files
            [filepath,name,ext] = fileparts(flname);
            if isempty(Filter)
                Files = dir([filepath '/' name '/*.tif']);
            else
                Files = dir([filepath '/' name '/*' Filter '*.tif']);
            end
            num_slices = numel(Files);
            if num_slices == 1
                Mode2D = 1;
                BlckSizeZ = 1;
            else
                Mode2D = 0;
            end
            if isempty(Files)
                defaultValues = {''};
                titleBar = 'No image matching filter in folder';
                userPrompt = {'Optional filename substring: '};
                caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                if ~isempty(caUserInput)
                    Filter = caUserInput{1};
                end
                flname = uigetdir('Images/','Select a 3D TIFF image folder');
                [filepath,name,ext] = fileparts(flname);
                aname = [filepath '\' name '_ann.tif'];
            else
                select = 0;
            end
        end
        
        info = imfinfo([flname '/' Files(1).name]);
        OrigWidth = info(1).Width;
        OrigHeight = info(1).Height;
        
        %% Downsampling?
        if OrigWidth*OrigHeight*num_slices > SizeLimit
            RecDownSamp = ceil(sqrt(OrigWidth*OrigHeight*num_slices/SizeLimit));
            str = inputdlg('Recommended downsampling factor','Large image detected',1,{num2str(RecDownSamp)});
            DownSamp = max([str2double(str{1}) 1]);
        else
            DownSamp = 1;    
        end
        
        %% Initialize arrays
        Width = floor(info(1).Width/DownSamp);
        Height = floor(info(1).Height/DownSamp);
        if info(1).BitDepth == 8
            I = uint8(zeros(Height,Width,num_slices));
        else
            I = uint16(zeros(Height,Width,num_slices));
        end
        
        %% Load images
        for kf = 1:num_slices
            h = waitbar(kf/num_slices);
            slice = imread([flname '/' Files(kf).name], 'PixelRegion',{[1+floor(DownSamp/2),DownSamp,Height*DownSamp],[1+floor(DownSamp/2),DownSamp,Width*DownSamp]});
            %% Color mode
            if size(slice,3) > 1
                if kf == 1
                    disp('Color images not supported, loading RGB average');
                end
                slice = round(mean(slice,3));
            end
            I(:,:,kf) = slice;
        end  
    end
    if ishandle(h)
        close(h);
    end

    %% Create empty annotation mask
    Lbl = uint8(zeros(Height,Width,num_slices));
    
    %% Current default annotation name
    if FolderMode == 1
        [filepath,name,ext] = fileparts(fname);
        aname = [filepath '\' name '.png'];
        anameref = [filepath '\' name '_ref.png'];
    end  
   
    %% Pad I and L to a block size multiple 
    disp(['Image size: ' mat2str(size(I))]);
    if (BlckSizeXY>1 || BlckSizeZ>1 || sum(Scales)>1) && (FolderMode == 0)
        disp('Padding image to a size compatible with pyramid classifier...');
    end
    if Mode2D == 0
        PadY = mod(BlckSizeXY*Scales(end)-mod(size(I,1),BlckSizeXY*Scales(end)),BlckSizeXY*Scales(end));
        PadX = mod(BlckSizeXY*Scales(end)-mod(size(I,2),BlckSizeXY*Scales(end)),BlckSizeXY*Scales(end));
        PadZ = mod(BlckSizeZ*Scales(end)-mod(size(I,3),BlckSizeZ*Scales(end)),BlckSizeZ*Scales(end));
        I = padarray(I,[PadY PadX PadZ],'post');
        Lbl = padarray(Lbl,[PadY PadX PadZ],'post'); 
        disp(['Size after padding: ' mat2str(size(I))]);
    else
        PadY = mod(BlckSizeXY*Scales(end)-mod(size(I,1),BlckSizeXY*Scales(end)),BlckSizeXY*Scales(end));
        PadX = mod(BlckSizeXY*Scales(end)-mod(size(I,2),BlckSizeXY*Scales(end)),BlckSizeXY*Scales(end));
        I = padarray(I,[PadY PadX],'post');
        Lbl = padarray(Lbl,[PadY PadX],'post'); 
        disp(['Size after padding: ' mat2str(size(I))]);
    end
     
    %% Display image
    imviewer = imtool3DLbl_mod(I,[0 0 1 1],handle);
    setCurrentSlice(imviewer,1+floor(num_slices/2));
    setMask(imviewer,Lbl);
    clear Lbl;
    Annotate = true;
    Interp = false;

    %% Main loop
    while Annotate
        try
            
        %% Show color coded annotations
        if ShowAnnot == true
            setAlpha(imviewer,Alpha);
        else
            setAlpha(imviewer,0);
        end
        
        %% Wait for command
        if FolderMode == 0
            set(handle,'Name','(P)aint   c(Y)linder   (W)and   (B)lob   s(K)eleton   (C)lassify   (D)L-batch      ||     (E)rase    (T)ransfer    (F)it    (G)rab    (S)ave   (L)oad   ||      (M)ask   (Z)proj   (V)olume      ||     (N)ew   (O)ptions   (H)elp   (Q)uit','NumberTitle','off');
            figure(handle);
        else
            set(handle,'Name','(P)aint   c(Y)linder   (W)and   (B)lob   s(K)eleton   (C)lassify   (D)L-batch      ||     (E)rase    (T)ransfer    (F)it    (G)rab    (S)ave   (L)oad   ||      (M)ask   (Z)proj   (V)olume      ||     ne(X)t   (N)ew   (O)ptions   (H)elp   (Q)uit)','NumberTitle','off');
            figure(handle);
        end
       
        %% Wait for user input
        KeyRead = '';
        while isempty(KeyRead)
            pause(0.05);
            Mode = KeyRead;
        end
        
        switch Mode

            %case 'i' %% Invert mask
            %    
            %    if Inv == 0
            %        imviewer = imtool3DLbl_mod(max(I(:))-I,[0 0 1 1],handle);
            %        Inv = 1;
            %    else
            %        imviewer = imtool3DLbl_mod(I,[0 0 1 1],handle);
            %        Inv = 0;
            %    end
                
            case 'n' %% New image (restart)
                
                answer = questdlg('Load new image?');
                if answer(1) == 'Y'
                    LOBSTER_Annotator;
                end
                
            case 'y' %% Cylinder
                
                old = imviewer.mask;
                stop = 0;
                Pos = [];
                set(handle,'Name','Cylinder: Set diameter','NumberTitle','off'); 
                h1 = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                position = h1.Position;
                h2 = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                position2 = h2.Position;
                if ~isempty(position) && ~isempty(position2)
                    delete(h1);
                    delete(h2);
                    x = position(1);
                    y = position(2);
                    x2 = position2(1);
                    y2 = position2(2);
                    Rad = ceil(sqrt((x2-x).^2+(y2-y).^2)/2);
                    disp(['Estimated radius: ' num2str(Rad)]);
                    disp('Lay axis nodes');
                    cnt = 0;
                    while stop == 0
                        set(handle,'Name','Lay nodes, ESC to complete','NumberTitle','off');
                        figure(handle);
                        cnt = cnt+1;
                        h(cnt) = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]); 
                        position = h(cnt).Position;
                        if ~isempty(position)
                            x = floor(position(1)+0.5);
                            y = floor(position(2)+0.5);
                            z = round(getCurrentSlice(imviewer));
                            Pos = [Pos;x y z];
                        else
                            stop = 1;
                        end
                    end
                    for i = 1:cnt
                        delete(h(i));
                    end
                    if size(Pos,1) >= 2
                        imviewer.mask = LayTube(imviewer.mask,Pos,Rad,LUT(3),ZRatio);
                    else
                        uiwait(msgbox('Minimum two nodes should be set!'));
                    end
                else
                    uiwait(msgbox('Two points should be set!'));
                end
                
            case 'v' %% 3D render
                
                LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle);
                
            case 'h'

                winopen('LOBSTER_Annotator.pdf');
                
            case 'g' %% Grab object
                
                old = imviewer.mask;
                set(handle,'Name','Click object to export','NumberTitle','off');
                h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                position = h.Position;
                delete(h);
                if ~isempty(position) 
                    x = floor(position(1)+0.5);
                    y = floor(position(2)+0.5);
                    z = round(getCurrentSlice(imviewer));
                    if old(y,x,z) > 0
                        M = uint8(regionGrow(old == old(y,x,z),y,x,z,0,1,Inf));
                        disp(['Estimated volume: ' num2str(sum(M(:)>0)) ' voxels']);
                        cp = bwconncomp(M>0);
                        stats = regionprops(cp,'BoundingBox');
                        BB = ceil(stats.BoundingBox);
                        expdir = uigetdir('Images\','Select an empty folder to export the object region');
                        h = waitbar(0,'Exporting object');
                        if ImageSeries == 1
                            defaultValues = {''};
                            titleBar = 'Filename filter substitution (leave empty for none)';
                            userPrompt = {'Optional substitution string: '};
                            caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                            if ~isempty(caUserInput)
                                Filter2 = caUserInput{1};
                            end
                        end
                        for kf = BB(3) : BB(3)+BB(6)-1
                            h = waitbar((kf-BB(3))/BB(6));
                            if ImageSeries == 0
                                slice = imread(fname, kf, 'PixelRegion',{[1+floor((BB(2)-1)*DownSamp),1,1+floor((BB(2)+BB(5)-1)*DownSamp)],[1+floor((BB(1)-1)*DownSamp),1,1+floor((BB(1)+BB(4)-1)*DownSamp)]});
                            else
                                rname = [flname '\' Files(kf).name];
                                if ~isempty(Filter2)
                                    rname = strrep(rname,Filter,Filter2);   
                                end
                                slice = imread(rname, 'PixelRegion',{[1+floor((BB(2)-1)*DownSamp),1,1+floor((BB(2)+BB(5)-1)*DownSamp)],[1+floor((BB(1)-1)*DownSamp),1,1+floor((BB(1)+BB(4)-1)*DownSamp)]});
                            end
                            if isa(slice, 'uint16')
                                msk = uint16(imviewer.mask(BB(2):BB(2)+BB(5),BB(1):BB(1)+BB(4),kf)>0);
                            else
                                msk = uint8(imviewer.mask(BB(2):BB(2)+BB(5),BB(1):BB(1)+BB(4),kf)>0);
                            end
                            slice = slice.*(imresize(msk,size(slice),'nearest'));
                            imwrite(slice, [expdir '\Export_Z' sprintf('%04d', kf) '.tif']);
                        end
                        if ishandle(h)
                            close(h);
                        end
                        clear M;
                    else
                        disp('Seed voxel is background');
                    end
                end
                
            case 't' %% Transfer annotations
                
                defaultValues = {'1','2','3','4'};
                titleBar = 'Transfer';
                userPrompt = {'Class 1 to class ','Class 2 to class','Class 3 to class','Class 4 to class'};
                caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                old = imviewer.mask;
                if ~isempty(caUserInput)
                    c1 = str2double(caUserInput{1});
                    c2 = str2double(caUserInput{2});
                    c3 = str2double(caUserInput{3});
                    c4 = str2double(caUserInput{4});
                    Lbl = imviewer.mask;
                    if c1+1 ~= 2
                        Lbl(imviewer.mask==LUT(2)) = LUT(c1+1);
                    end
                    if c2+1 ~= 3
                        Lbl(imviewer.mask==LUT(3)) = LUT(c2+1);
                    end
                    if c3+1 ~= 4
                        Lbl(imviewer.mask==LUT(4)) = LUT(c3+1);
                    end
                    if c4+1 ~= 5
                        Lbl(imviewer.mask==LUT(5)) = LUT(c4+1);
                    end
                    setMask(imviewer,Lbl);
                    disp('Done!');
                    clear Lbl;
                end
                
            case 'o' %% Options
                
                titleBar = 'Set options';
                if FolderMode == 0
                    defaultValues = {num2str(ZRatio),num2str(RunProj),num2str(ZSpan),num2str(PointSize),num2str(Alpha)};
                    userPrompt = {'ZRatio (0.25 - 2.5)','Z projection depth (3-127 slices)','3D render depth (8-Inf)','3D render point size (1-5)','Mask opacity (0.1-1)'};
                else
                    defaultValues = {num2str(ZRatio),num2str(RunProj),num2str(ZSpan),num2str(PointSize),num2str(Alpha),num2str(AnnSource)};
                    userPrompt = {'ZRatio (0.25 - 2.5)','Z projection depth (3-127 slices)','3D render depth (8-Inf)','3D render point size (1-5)','Mask opacity (0.1-1)','Load Reference annotations? (0/1)'};
                end
                caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                if ~isempty(caUserInput)
                    ZRatio = boundvar(str2double(caUserInput{1}),0.25,2.5,1);
                    RunProj = boundvar(str2double(caUserInput{2}),3,127,0);
                    ZSpan = boundvar(str2double(caUserInput{3}),8,Inf,Inf);
                    PointSize = boundvar(str2double(caUserInput{4}),1,5,1);
                    Alpha = boundvar(str2double(caUserInput{5}),0.1,1,1);
                    if FolderMode == 1
                        AnnSource = boundvar(str2double(caUserInput{6}),0,1,0);
                    end
                end
                if RunProjViewer > 0
                    RunProjViewer = RunProj;
                end
                
            case 'x' %% Next image (folder mode)
                
                if FolderMode == 1
                    Annotate = false;
                end
                
            case 'q' %% Exit
                
                answer = questdlg('Confirm exit?');
                if answer(1) == 'Y'
                    close all;
                    error('Stopped by user');
                end
                
            case 'u' %% Undo
                
                if exist('old','var')
                    setMask(imviewer,uint8(old));
                    clear old;
                end
                
            case 'z' %% Toggle local Z projection on/off
                
                projmode = ~projmode;
                if projmode == 1
                    RunProjViewer = RunProj;
                else
                    RunProjViewer = 0;
                end      
                   
            case 'e' %% Erase annotations
                
                answer = questdlg('Erase annotations?');
                if answer(1) == 'Y'
                    %% Clear annotations
                    old = imviewer.mask;
                    setMask(imviewer,imviewer.mask*0);
                end
                
            case 'm' %% Toggle annotations on/off
                
                ShowAnnot = ~ShowAnnot;
                if ShowAnnot == false
                    setAlpha(imviewer,0);
                else
                    setAlpha(imviewer,Alpha);
                end
                
            case 'd' %% Deep learning batch processing (2D images only)
                
                answer = questdlg('Batch process image folder (load results with "l")?');
                
                if answer(1) == 'Y' 
                
                    disp('Batch process image folder...');
                    if FolderMode == 1
                        h = waitbar(0,'Batch process image folder...');
                        UNET_pixclass(selpath,[Height Width],LUT(2),LUT(3));
                        if ishandle(h)
                            close(h);
                        end
                        disp('Done!');
                        [filepath,name,ext] = fileparts(fname);
                        aname = [filepath '\' name '.png'];
                    else
                        disp('This mode is only compatible with folder of 2D TIFF images');
                    end
                
                end
                
            case 's' %% Save annotations
                
                disp('Saving annotations...');
                Lbl = imviewer.mask(1:Height,1:Width,1:num_slices);
                Save = 1;
                if exist(aname,'file')
                    answer = questdlg('File already exists, overwrite?');
                    if answer(1) == 'Y'
                        delete(aname);
                    else
                        Save = 0;
                    end
                end
                if Save == 1
                    h = waitbar(0,'Saving annotations');
                    if DownSamp > 1
                        imwrite(imresize(uint8(Lbl(:,:,1)),[OrigHeight OrigWidth],'nearest'), aname);
                    else
                        imwrite(uint8(Lbl(:,:,1)), aname);
                    end
                    for kf = 2:num_slices
                        h = waitbar(kf/num_slices);
                        if DownSamp > 1
                            imwrite(imresize(uint8(Lbl(:,:,kf)),[OrigHeight OrigWidth],'nearest'), aname, 'WriteMode', 'append');
                        else
                            imwrite(uint8(Lbl(:,:,kf)), aname, 'WriteMode', 'append');
                        end
                    end
                    if ishandle(h)
                        close(h);
                    end
                    clear Lbl;
                end
                
            case 'l' %% Load annotations
                
                if ~exist(aname,'file')
                    [file,path] = uigetfile('*.tif','Select a 2D/3D annotation TIFF file');
                    aname_load = [path file];
                else
                    aname_load = aname;
                    if FolderMode == 1 && AnnSource == 1
                        aname_load = anameref;
                    end
                    file = 1;
                end
                if file ~= 0
                    disp('Loading annotations...');
                    h = waitbar(0,'Loading annotations...');
                    info = imfinfo(aname_load);
                    Width_ann = floor(info(1).Width/DownSamp);
                    Height_ann = floor(info(1).Height/DownSamp);
                    num_slices_ann = numel(info);
                    Lbl = uint8(zeros(size(imviewer.mask)));
                    if num_slices_ann <= num_slices && Width_ann <= size(imviewer.mask,2) && Height_ann <= size(imviewer.mask,1)
                        for kf = 1:num_slices_ann
                            waitbar(kf/num_slices);
                            if DownSamp == 1
                                if num_slices_ann > 1
                                    slice = imread(aname_load, kf);
                                else
                                    slice = imread(aname_load);
                                end
                            else
                                if num_slices_ann > 1
                                    slice = imread(aname, kf, 'PixelRegion',{[1,DownSamp,Height*DownSamp],[1,DownSamp,Width*DownSamp]});
                                else
                                    slice = imread(aname, 'PixelRegion',{[1,DownSamp,Height*DownSamp],[1,DownSamp,Width*DownSamp]});
                                end
                            end
                            Lbl(1:Height_ann,1:Width_ann,kf) = slice;
                        end
                        if FolderMode == 1
                            if max(Lbl(:)) == 255
                                Lbl(Lbl==0) = LUT(2);
                                Lbl(Lbl==255) = LUT(3);
                            end
                        end
                        setMask(imviewer,Lbl);
                    else
                        uiwait(msgbox('Incompatible annotation size!'));
                    end
                    if ishandle(h)
                        close(h);
                    end
                    disp('Done!');
                else
                    disp('No annotation file selected');
                end
                
            case 'p' %% Paint
                
                Paint = true;
                while Paint == true 
                    
                    if Polarity == 1
                        set(handle,'Name','Paint: (0-3)   (I)nterp_Class1   (B)rush:bright           ||           (M)ask   (Z)proj   (V)olume    (U)ndo    e(X)it','NumberTitle','off');
                        figure(handle);
                    else
                        set(handle,'Name','Paint: (0-3)   (I)nterp_Class1   (B)rush:dark             ||           (M)ask   (Z)proj   (V)olume    (U)ndo    e(X)it','NumberTitle','off');
                        figure(handle);
                    end
                   
                    %% Show color coded annotations
                    if ShowAnnot == true
                        setAlpha(imviewer,Alpha);
                    else
                        setAlpha(imviewer,0);
                    end
                    
                    %% Wait for user input
                    KeyRead = '';
                    while isempty(KeyRead)
                        pause(0.05);
                        currkey = KeyRead;
                    end
                    
                    switch currkey
                        
                        case 'x' %% Exit
                            
                            Paint = false;
                            
                        case 'm' %% Toggle mask on/off
                            
                            ShowAnnot = ~ShowAnnot;
                            if ShowAnnot == false
                                setAlpha(imviewer,0);
                            else
                                setAlpha(imviewer,Alpha);
                            end
                            
                        case 'b' %% Toggle magic brush polarity
                            
                            Polarity = ~Polarity;
                            
                        case 'i' %% Interpolate slices from class 1 to class 2 objects
                            
                            old = imviewer.mask;
                            disp('Interpolating');
                            h = waitbar(1,'Interpolating');
                            if size(old,3)>1
                                sliceindx = [];
                                for i = 1:size(old,3)
                                    currentslice = (old(:,:,i)==LUT(2));
                                    if sum(currentslice(:)>0)
                                        sliceindx = [sliceindx i]; 
                                    end
                                end
                                for i = 1:numel(sliceindx)-1
                                    z1 = sliceindx(i);
                                    z2 = sliceindx(i+1);
                                    vol = uint8(interp_shape(single(imfill(old(:,:,z2)==LUT(2),'holes')),single(imfill(old(:,:,z1)==LUT(2),'holes')),1+z2-z1))*LUT(3);
                                    buf = imviewer.mask(:,:,z1:z2);
                                    buf(vol>0) = 0;
                                    buf = buf+vol;
                                    buf(buf==LUT(2)) = 0;
                                    imviewer.mask(:,:,z1:z2) = buf;
                                end
                            else
                                vol = uint8(old == LUT(2));
                                vol = uint8(imfill(vol,'holes')*LUT(3));
                                buf = old;
                                buf(vol>0) = 0;
                                buf = buf+vol;
                                imviewer.mask = buf;
                            end
                            clear buf;clear vol;
                            notify(imviewer,'maskChanged');
                            disp('Done!');
                            if ishandle(h)
                                close(h);
                            end
                            
                        case 'z' %% Toggle local Z projection
                            projmode = ~projmode;
                            if projmode == 1
                                RunProjViewer = RunProj;
                            else
                                RunProjViewer = 0;
                            end
                        
                        case 'v' %% 3D render
                
                            LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle);    
                            
                        case 'u' %% Undo
                            
                            if exist('old','var')
                                setMask(imviewer,uint8(old));
                                clear old;
                            end
                            
                        otherwise %% Paint Class 1-3
                            
                            old = imviewer.mask;
                            c = str2double(currkey);
                            if ~isempty(c)
                                if c>=0 && c<= 3
                                    if(c>0)
                                        LastColor = LUT(1+c);
                                    end
                                    set(handle,'Name',['Class' num2str(c)],'NumberTitle','off');
                                    h = drawfreehand('DrawingArea',[0 0 size(I,1) size(I,2)]);
                                    coords = round(h.Position);
                                    delete(h);
                                    z = round(getCurrentSlice(imviewer));
                                    buf = uint8(imviewer.mask(:,:,z));
                                    if c>0
                                        slice = uint8(poly2mask(coords(:,1),coords(:,2),size(I,1),size(I,2)))*LUT(c+1);
                                        buf(slice>0) = 0;
                                        imviewer.mask(:,:,z) = buf+slice;
                                    else
                                        slice = uint8(poly2mask(coords(:,1),coords(:,2),size(I,1),size(I,2)));
                                        buf(slice>0) = 0;
                                        imviewer.mask(:,:,z) = buf;
                                    end
                                    notify(imviewer,'maskChanged');
                                    clear buf;
                                end
                            end
                        end
                end
                if exist('slice2','var')
                    clear slice2;
                end       
            
            case 'w' %% Wand
                  
                Wand = true;
                while Wand == true
                    
                    set(handle,'Name','Wand: (0-3)   (C)ut   (S)plit_Class1   (J)oin_Class2   (R)emove   (F)ix   (O)ptions          ||          (M)ask   (Z)proj   (V)olume    (U)ndo   e(X)it','NumberTitle','off');
                    figure(handle);
                    
                    %% Show color coded annotations
                    if ShowAnnot == true
                        setAlpha(imviewer,Alpha);
                    else
                        setAlpha(imviewer,0);
                    end
                    
                    %% Wait for user input
                    KeyRead = '';
                    while isempty(KeyRead)
                        pause(0.05);
                        currkey = KeyRead;
                    end
                    
                    switch currkey
                        
                        case 'x'  %% Exit
                            Wand = false;
                            
                        case 'z' %% Local Z projection toggle
                            projmode = ~projmode;
                            if projmode == 1
                                RunProjViewer = RunProj;
                            else
                                RunProjViewer = 0;
                            end
                        
                        case 'v' %% 3D render
                
                            LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle);    
                            
                        case 'u' %% Undo
                            if exist('old','var')
                                setMask(imviewer,uint8(old));
                                clear old;
                            end
                            
                        case 'o' %% Options
                            
                            defaultValues = {num2str(wandblur), num2str(wandtol), num2str(wanddst), num2str(wanddim), num2str(wandmode), num2str(splitmode), num2str(gaussraddstmap), num2str(Min2DHolesArea)};
                            titleBar = 'Set parameters for wand / split';
                            userPrompt = {'Blur (0-7)', 'Tolerance (0-1)', 'Wand distance (5-Inf)', 'Wand dim (0:2D, 1:3D)', 'Threshold (0:low/up, 1:low)','Splitting: algorithm (0: regiondescent, 1: watershed)', 'Splitting: Distance map blur (0-15)','Kept 2D holes minimum area (0-1000)'};
                            caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                            if ~isempty(caUserInput)
                                wandblur = boundvar(str2double(caUserInput{1}),0,7,1);
                                wandtol = boundvar(str2double(caUserInput{2}),0,1,0.5);
                                wanddst = boundvar(str2double(caUserInput{3}),5,Inf,5);
                                wanddim = boundvar(str2double(caUserInput{4}),0,1,1)>0;
                                wandmode = boundvar(str2double(caUserInput{5}),0,1,0)>0;
                                splitmode = boundvar(str2double(caUserInput{6}),0,1,1);
                                gaussraddstmap = boundvar(str2double(caUserInput{7}),0,15,2);
                                Min2DHolesArea = boundvar(str2double(caUserInput{8}),0,1000,15);
                                %growdst = boundvar(str2double(caUserInput{8}),1,100,100);
                            end
                            
                        case 'r' %% Remove object
                            
                            old = imviewer.mask;
                            set(handle,'Name','Click object to remove','NumberTitle','off');
                            h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                            position = h.Position;
                            delete(h);
                            if ~isempty(position) 
                                x = floor(position(1)+0.5);
                                y = floor(position(2)+0.5);
                                z = round(getCurrentSlice(imviewer));
                                if old(y,x,z) > 0
                                    M = uint8(regionGrow(old == old(y,x,z),y,x,z,0,1,Inf));
                                    setMask(imviewer,imviewer.mask-uint8(old(y,x,z)*M));
                                    clear M;
                                else
                                    disp('Seed voxel is background');
                                end
                            end
                            
                        case 'f' %% Keep only selected object
                            
                            old = imviewer.mask;
                            set(handle,'Name','Click object to keep','NumberTitle','off');
                            h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                            position = h.Position;
                            delete(h);
                            if ~isempty(position) 
                                x = floor(position(1)+0.5);
                                y = floor(position(2)+0.5);
                                z = round(getCurrentSlice(imviewer));
                                if old(y,x,z) > 0
                                    M = uint8(regionGrow(old == old(y,x,z),y,x,z,0,1,Inf));
                                    setMask(imviewer,uint8(old(y,x,z)*M));
                                    clear M;
                                else
                                    disp('Seed voxel is background');
                                end
                            end
                            
                        case 'c' %% Cut
                            
                             old = imviewer.mask;
                             disp('Draw cut line in top slice');
                             set(handle,'Name','Draw cut line in top slice (hold left, release)','NumberTitle','off');
                             h = drawline('DrawingArea',[0 0 size(I,1) size(I,2)]);
                             position = h.Position;
                             delete(h);
                             if ~isempty(position)
                                 [x y] = bresenham(position(1,1),position(1,2),position(2,1),position(2,2));
                                 inds = y+x*size(old,1);
                                 inds = [inds inds+1 inds+size(old,1)];
                                 s1 = getCurrentSlice(imviewer);
                                 if Mode2D == 0
                                    disp('Mark point in bottom slice');
                                    set(handle,'Name','Mark a point in bottom slice','NumberTitle','off');
                                    h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                                    position = h.Position;
                                    delete(h);  
                                 end
                                 if ~isempty(position)
                                     s2 = getCurrentSlice(imviewer);
                                     smin = min([s1 s2]);
                                     smax = max([s1 s2]);
                                     for i = smin:smax
                                        imviewer.mask(inds+(i-1)*size(old,1)*size(old,2)) = 0;
                                     end
                                 end
                             end
                             
                        case 's' %% Split Class 1 objects with watershed or region descent
                            
                            disp('Splitting Class 1 objects...');
                            old = imviewer.mask;
                            vol = (old==LUT(2));
                            if Min2DHolesArea>0
                                for i = 1:size(vol,3)
                                    vol(:,:,i) = vol(:,:,i)|(~(~vol(:,:,i) & bwareaopen(~vol(:,:,i),Min2DHolesArea)));
                                end
                            end
                            D = bwdist(~vol);
							D = imgaussfilt(D,gaussraddstmap);
                            if splitmode == 0
                                disp('region descent...');
                                h = waitbar(1,'Splitting with region descent');
                                S = (imdilate(D,ones(3,3,3))==D);
                                S(~vol) = 0;
                                D(1,:,:) = Inf;D(end,:,:) = Inf;D(:,1,:) = Inf;D(:,end,:) = Inf;D(:,:,1) = Inf;D(:,:,end) = Inf;
                                S(1,:,:) = 0;S(end,:,:) = 0;S(:,1,:) = 0;S(:,end,:) = 0;S(:,:,1) = 0;S(:,:,end) = 0;
                                dlt = size(S,1)*size(S,2);
                                offs = [-1 1 -size(S,1) size(S,1) -1-size(S,1) 1-size(S,1) -1+size(S,1) 1+size(S,1)... 
                                        -dlt-1 -dlt+1 -dlt-size(S,1) -dlt+size(S,1) -dlt-1-size(S,1) -dlt+1-size(S,1) -dlt-1+size(S,1) -dlt+1+size(S,1)...
                                        dlt-1 dlt+1 dlt-size(S,1) dlt+size(S,1) dlt-1-size(S,1) dlt+1-size(S,1) dlt-1+size(S,1) dlt+1+size(S,1)].';
                                Seeds = single(find(S>0));
                                Regs = num2cell(Seeds);
                                L = single(zeros(size(S)));
                                L(Seeds) = 1:numel(Seeds);
                                for it = 1:splitdst
                                    for r = 1:numel(Regs)
                                        LastDs = D(Regs{r});
                                        Cands = repmat(Regs{r},numel(offs),1)+repmat(offs,1,numel(Regs{r}));
                                        ValCands = (D(Cands)<repmat(LastDs,numel(offs),1))&(D(Cands)>1)&(L(Cands)==0);
                                        Regs{r} = single(Cands(ValCands));
                                        L(Cands(ValCands)) = r;
                                        Regs{r} = unique((Regs{r}(:))).';
                                    end
                                    if sum(cellfun(@isempty,Regs)) == numel(Regs)
                                        it = growdst;
                                    end
                                end
                                Edg = imdilate(L,ones(3,3,3))-L;
                                L(Edg>0) = 0;
                            else
                                disp('watershed...');
                                h = waitbar(1,'Splitting with watershed');
                                D(~vol) = Inf;
                                L = watershed(-D);
                                L(~vol) = 0;
                            end
                            clear D;
                            buf = old;
                            vol = uint8(single(LUT(3))*(L>0));
                            clear L;
                            buf(vol>0) = 0;
                            buf = buf+vol;
                            buf(buf==LUT(2)) = 0;
                            imviewer.mask = buf;
                            clear buf;clear vol;
                            if ishandle(h)
                                close(h);
                            end
                            disp('Done!');
                            
                        case 'j' %% Join objects
                            
                            old = imviewer.mask;
                            set(handle,'Name','Mark first object','NumberTitle','off');
                            h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                            position = h.Position;
                            delete(h);
                            if ~isempty(position)
                                x = floor(position(1)+0.5);
                                y = floor(position(2)+0.5);
                                z = round(getCurrentSlice(imviewer));
                                if old(y,x,z) == LUT(3)
                                    M1 = uint8(regionGrow(old == LUT(3),y,x,z,0,1,Inf));
                                else
                                    disp('Invalid first object');
                                end
                                set(handle,'Name','Mark second object','NumberTitle','off');
                                h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                                position2 = h.Position;
                                delete(h);
                                x2 = floor(position2(1)+0.5);
                                y2 = floor(position2(2)+0.5);
                                z2 = round(getCurrentSlice(imviewer));
                                if old(y2,x2,z2) == LUT(3)
                                    M2 = uint8(regionGrow(old == LUT(3),y2,x2,z2,0,1,Inf));
                                else
                                    disp('Invalid second object');
                                end
                                if exist('M1','var') && exist('M2','var')
                                    M1 = M1|M2;
                                    M1 = imclose(M1,ones(3,3,3));
                                    setMask(imviewer,uint8(single(LUT(3))*(imviewer.mask|M1)));
                                end
                                clear M1;clear M2;
                            end
                            
                        case 'm' %% Toggle mask on/off
                            
                            ShowAnnot = ~ShowAnnot;
                            if ShowAnnot == false
                                setAlpha(imviewer,0);
                            else
                                setAlpha(imviewer,Alpha);
                            end
                            
                        otherwise %% Wand and store to Class 1-3
                            
                            c = str2double(currkey);
                            if ~isempty(c)
                                if c>=0 && c<= 3
                                    set(handle,'Name',['Class' num2str(c)],'NumberTitle','off');
                                    h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                                    position = h.Position;
                                    delete(h);
                                    if ~isempty(position)
                                        x = floor(position(1)+0.5);
                                        y = floor(position(2)+0.5);
                                        z = round(getCurrentSlice(imviewer));
                                        if wanddim == 0
                                            if wandblur>0
                                                M = regionGrow2D(imgaussfilt(I(:,:,z),wandblur),y,x,wandtol,wandmode,wanddst);
                                            else
                                                M = regionGrow2D(I,y,x,wandtol,wandmode,wanddst);
                                            end
                                            old = imviewer.mask;
                                            buf = imviewer.mask(:,:,z);
                                            imviewer.mask(:,:,z) = buf.*uint8(M==0)+uint8(M>0)*LUT(1+c);
                                            notify(imviewer,'maskChanged');
                                        else
                                            h = waitbar(0,'Exploring from seed');
                                            if wandblur>0
                                                M = regionGrow(imgaussfilt(I,wandblur),y,x,z,wandtol,wandmode,wanddst);
                                            else
                                                M = regionGrow(I,y,x,z,wandtol,wandmode,wanddst);
                                            end
                                            if ishandle(h)
                                                close(h);
                                            end
                                            old = imviewer.mask;
                                            imviewer.mask = old.*uint8(M==0)+uint8(M>0)*LUT(1+c);
                                        end
                                    end
                                end
                            end
                    end
                end
                
            case 'b' %% Blob 
                  
                 %% Main loop
                 Seed = true;
                 while Seed == true
                     
                    set(handle,'Name','Blob: (S)eed   (A)dd   (R)emove   (F)ill   (O)ptions           ||           (M)ask   (Z)proj   (V)olume    (U)ndo   e(X)it','NumberTitle','off');
                    figure(handle);
                    
                    %% Show color coded annotations
                    if ShowAnnot == true
                        setAlpha(imviewer,Alpha);
                    else
                        setAlpha(imviewer,0);
                    end  
                    
                    %% Wait for user input
                    KeyRead = '';
                    while isempty(KeyRead)
                        pause(0.05);
                        currkey = KeyRead;
                    end
                    
                    switch currkey
                        
                        case 'x' %% Exit
                            
                            Seed = false;
                            
                        case 'o' %% Options
                            
                            defaultValues = {num2str(Sig), num2str(LocalMaxBox), num2str(MinLoG), num2str(FillMode)};
                            titleBar = 'Set parameters for blob seeding';
                            userPrompt = {'Sigma (0.5-8)', 'LocalMaxBox (3-25)', 'MinLog (0.01-1)', 'FillMode (0:watershed 1:propagate)'};
                            caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                            if ~isempty(caUserInput)
                                Sig = boundvar(str2double(caUserInput{1}),0.5,8,4);
                                LocalMaxBox = boundvar(str2double(caUserInput{2}),3,25,7);
                                MinLoG = boundvar(str2double(caUserInput{3}),0.01,1,0.1);    
                                FillMode = boundvar(str2double(caUserInput{4}),0,1,0);
                            end
                            
                        case 's' %% Seed
                            
                             disp('Automatic seeding...');
                             h = waitbar(1,'Automatic seeding');
                             old = imviewer.mask; 
                             if Mode2D == 0
                                params.Sigmas = Sig*[1 1 1];
                                params.LocalMaxBox = LocalMaxBox*[1 1 1];
                                params.MinLoG = MinLoG;
                                S = fxg_sLoG3DLocMax3D(I, params);
                                S = imdilate(S,ones(3,3,3));
                             else
                                params.GRad = Sig;
                                params.LocalMaxBox = LocalMaxBox*[1 1];
                                params.MinLoG = MinLoG;
                                S = fxg_sLoGLocMax(I, params);
                                S = imdilate(S,ones(3,3));
                             end
                             setMask(imviewer,uint8(single(LUT(3))*(S>0)));
                             clear S;clear buf;
                             if ishandle(h)
                                 close(h);
                             end
                             disp('Done!');
                             
                        case 'f' %% Fill from seeds
                            
                             old = imviewer.mask;
                             S = (imviewer.mask == LUT(3));
                             if sum(S(:))>0
                                  if FillMode == 1  
                                      disp('Growing seeds with propagate...');
                                      h = waitbar(1,'Growing seeds with propagate');
                                      params.Power = 1;
                                      params.AnalyzeCC = 1;
                                      if Mode2D == 0
                                            W = fxgs_lSeededPropagate3D_mod(single(I),S,true(size(S)),params);
                                      else
                                            params.Thr = 0;
                                            W = fxgs_lSeededPropagate(single(I),S,params);
                                      end
                                      W = uint8(single(LUT(3))*(imdilate(W,ones(3,3,3)) == W));
                                  else
                                      disp('Growing seeds with watershed...');
                                      h = waitbar(1,'Growing seeds with watershed');
                                      if Mode2D == 0
                                        If = single(imgaussfilt3(I,1));
                                        [Ifx,Ify,Ifz] = gradient(If);
                                        clear If;
                                        GradMag = sqrt(Ifx.^2+Ify.^2+Ifz.^2);
                                        clear Ifx;clear Ify;clear Ifz;
                                        GradMag = imimposemin(GradMag,(S>0),26);
                                      else
                                        If = single(imgaussfilt(I,1));
                                        [Ifx,Ify] = gradient(If);
                                        clear If;
                                        GradMag = sqrt(Ifx.^2+Ify.^2);
                                        clear Ifx;clear Ify;
                                        GradMag = imimposemin(GradMag,(S>0),8);
                                      end
                                      W = watershed(GradMag);
                                      clear GradMag;
                                      W = uint8(single(LUT(3))*(W>0));
                                  end
                                  W(1,:,:) = 0;W(end,:,:) = 0;W(:,1,:) = 0;W(:,end,:) = 0;
                                  setMask(imviewer,W);
                                  clear W;
                                  if ishandle(h)
                                    close(h);
                                  end
                                  disp('Done!');   
                             end
                             RunProjViewer = 0;
                             
                        case 'r' %% Remove object
                            
                            old = imviewer.mask;
                            set(handle,'Name','Mark object','NumberTitle','off');
                            h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                            position = h.Position;
                            delete(h);
                            if ~isempty(position) 
                                x = floor(position(1)+0.5);
                                y = floor(position(2)+0.5);
                                z = round(getCurrentSlice(imviewer));
                                if old(y,x,z) == LUT(3)
                                    M = uint8(regionGrow(old == LUT(3),y,x,z,0,1,Inf));
                                    setMask(imviewer,imviewer.mask-uint8(single(LUT(3)*M)));
                                    clear M;
                                else
                                    disp('Seed voxel is not from class 2');
                                end
                            end
                            
                        case 'a' %% Add seed
                            
                            old = imviewer.mask;
                            set(handle,'Name','Mark seed','NumberTitle','off');
                            h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                            position = h.Position;
                            delete(h);
                            if ~isempty(position)
                                x = floor(position(1)+0.5);
                                y = floor(position(2)+0.5);
                                z = round(getCurrentSlice(imviewer));
                                M = imviewer.mask;
                                M(y-1:y+1,x-1:x+1,z) = LUT(3);
                                setMask(imviewer,M);
                                clear M;
                            end
                            
                        case 'm' %% Toggle mask on/off
                            ShowAnnot = ~ShowAnnot;
                            if ShowAnnot == false
                                setAlpha(imviewer,0);
                            else
                                setAlpha(imviewer,Alpha);
                            end
                            
                        case 'z' %% Toggle Z projection on/off
                            projmode = ~projmode;
                            if projmode == 1
                                RunProjViewer = RunProj;
                            else
                                RunProjViewer = 0;
                            end
                        
                        case 'v' %% 3D render
                
                            LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle);    
                            
                        case 'u' %% Undo
                            
                            if exist('old','var')
                                setMask(imviewer,uint8(old));
                                clear old;
                            end
                        end
                 end
                    
            case 'k' %% Skeleton
                  
                 DmapComputed = 0;
                 clear Path;
                 Skeleton = true;
                 while Skeleton == true
                     
                    set(handle,'Name','Skeleton:   s(K)el_Class1   (T)race_Class1   (R)emove   (E)rase skeleton         ||              (M)ask   (Z)proj   (V)olume   (U)ndo   (O)ptions    e(X)it)','NumberTitle','off');
                    figure(handle);
                    
                    %% Show color coded annotations
                    if ShowAnnot == true
                        setAlpha(imviewer,Alpha);
                    else
                        setAlpha(imviewer,0);
                    end
                    
                    %% Wait for user input
                    KeyRead = '';
                    while isempty(KeyRead)
                        pause(0.05);
                        currkey = KeyRead;
                    end
                    
                    switch currkey
                        
                        case 'x' %% Exit
                            
                            Skeleton = false;
                            
                        case 'o' %% Options
                            
                            defaultValues = {num2str(PreCloseRad), num2str(Min2DHolesArea), num2str(MinVol), num2str(MinBrchLgth), num2str(PruneIter), num2str(Ignore4Way), num2str(brchmode)};
                            titleBar = 'Set parameters for skeletonization';
                            userPrompt = {'PreCloseRad (0-9)', '2D holes fill area (0-1000)', 'MinVol (0-1000)', 'Minimum branch length (0-64)', 'Pruning iterations (0-4)', 'Prune 4-way branching points', 'Trace mode (0:normal, 1:first point branch/end)'};
                            caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                            if ~isempty(caUserInput)
                                PreCloseRad = boundvar(str2double(caUserInput{1}),0,9,3);
                                Min2DHolesArea = boundvar(str2double(caUserInput{2}),0,1000,25);
                                MinVol = boundvar(str2double(caUserInput{3}),0,1000,50);
                                MinBrchLgth = boundvar(str2double(caUserInput{4}),0,64,7);
                                PruneIter = boundvar(str2double(caUserInput{5}),0,4,1);
                                Ignore4Way = boundvar(str2double(caUserInput{6}),0,1,0);
                                brchmode = boundvar(str2double(caUserInput{7}),0,1,0);
                            end
                            
                        case 'k' %% skeletonize
                            
                            if max(imviewer.mask(:))<LUT(6)
                                h = waitbar(1,'Operation in progress...');
                                disp('Skeletonizing class 1...');
                                old = imviewer.mask;
                                params.PreCloseRad = PreCloseRad;
                                params.Min2DHolesArea = Min2DHolesArea;
                                params.MinVol = MinVol;
                                new = fxm_kSkl3D(old == LUT(2), params);   
                                disp('Labeling skeleton...');
                                params.SklLbl = 1;
                                params.MinBrchLgth = MinBrchLgth;
                                params.MinBrchLgth2 = MinBrchLgth;
                                params.MaxIter = PruneIter;
                                params.Ignore4Way = Ignore4Way;
                                params.ZRatio = ZRatio;
                                [new] = fxk_kSklLbl3D(uint8(LUT(6))*uint8(new>0), params);
                                old = imviewer.mask;
                                buf = imviewer.mask;
                                buf(new>0) = new(new>0);
                                setMask(imviewer,uint8(buf));
                                clear buf;
                                disp('Done!');
                                if ishandle(h)
                                    close(h);
                                end
                            else
                                disp('Mask is already skeletonized, clear skeleton first!');
                            end
                            
                        case 'e' %% Clear skeleton
                            old = imviewer.mask;
                            new = old;
                            new(new>=LUT(5)) = LUT(2);
                            setMask(imviewer,new);
                            clear new;
                            
                        case 't' %% Trace branch manually
                            
                            if RunProjViewer == 0 
                                old = imviewer.mask;
                                SklClass = 1;
                                if DmapComputed ~= SklClass
                                    disp('Computing distance map...');
                                    h = waitbar(1,'Computing distance map');
                                    Dmap = single(zeros(size(imviewer.mask)));
                                    for i = 1:size(imviewer.mask,3)
                                        Dmap(:,:,i) = 1./bwdist((imviewer.mask(:,:,i)~=LUT(1+SklClass))&(imviewer.mask(:,:,i)~=LUT(6))&(imviewer.mask(:,:,i)~=LUT(7)));
                                    end
                                    disp('Done!');
                                    if ishandle(h)
                                        close(h);
                                    end
                                    DmapComputed = SklClass;
                                end
                                if brchmode == 1
                                    set(handle,'Name','Set first point','NumberTitle','off');
                                else
                                    set(handle,'Name','Set first point','NumberTitle','off');
                                end
                                h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                                position = h.Position;
                                delete(h);
                                if ~isempty(position)
                                    x1 = floor(position(1)+0.5);
                                    y1 = floor(position(2)+0.5);
                                    z1 = round(getCurrentSlice(imviewer));
                                    set(handle,'Name','Set second point','NumberTitle','off');
                                    h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                                    position2 = h.Position;
                                    delete(h);
                                    x2 = floor(position2(1)+0.5);
                                    y2 = floor(position2(2)+0.5);
                                    z2 = round(getCurrentSlice(imviewer));
                                    %mZ = max(min(z1,z2)-Tracedepth/2,1);
                                    %MZ = min(max(z1,z2)+Tracedepth/2,size(imviewer.mask,3));
                                    mZ = 1;
                                    MZ = size(I,3);
                                    FirstPoint = imviewer.mask(y1,x1,z1);
                                    SecondPoint = imviewer.mask(y2,x2,z2);
                                    if (FirstPoint>=LUT(2)) && (SecondPoint>=LUT(2))
                                        Gmap = graydist(Dmap(:,:,mZ:MZ),sub2ind(size(I),y2,x2,1+z2-mZ));
                                        Path = shortestpath(Gmap,y1,x1,1+z1-mZ,y2,x2,1+z2-mZ);
                                        imviewer.mask(Path+(mZ-1)*size(imviewer.mask,1)*size(imviewer.mask,2)) = LUT(6);
                                        if brchmode == 1
                                            if FirstPoint >= LUT(6)
                                                imviewer.mask(y1,x1,z1) = LUT(8);
                                            else
                                                imviewer.mask(y1,x1,z1) = LUT(7);
                                            end
                                            if SecondPoint >= LUT(6)
                                                imviewer.mask(y2,x2,z2) = LUT(8);
                                            else
                                                imviewer.mask(y2,x2,z2) = LUT(7);
                                            end  
                                        end
                                        notify(imviewer,'maskChanged');
                                    else
                                        disp('Both points should be non background');
                                    end
                                end
                            else
                                disp('Tracing not compatible with dynamic Z projection');
                            end
                            
                        case 'r' %% Remove branch manually
                            
                            old = imviewer.mask;
                            h = drawpoint('DrawingArea',[0 0 size(I,1) size(I,2)]);
                            position = h.Position;
                            delete(h);
                            if ~isempty(position)
                                x = floor(position(1)+0.5);
                                y = floor(position(2)+0.5);
                                z = round(getCurrentSlice(imviewer));
                                M = imviewer.mask;
                                if RunProjViewer > 0
                                    Zscan = M(y,x,z-round(RunProjViewer/2):z+round(RunProjViewer/2));
                                    cands = find(Zscan>=LUT(6))-round(numel(Zscan)/2);
                                    [minval,minind] = min(abs(cands));
                                    z = z + cands(minind);
                                end
                                if M(y,x,z) >= LUT(2)
                                    if M(y,x,z) == LUT(8)
                                        %% Click on branch point --> remove branch point
                                        M = (M == LUT(8));
                                        M = regionGrow26(M,y,x,z,0,Inf);
                                        setMask(imviewer,imviewer.mask-uint8(LUT(8)-LUT(6))*uint8(M));
                                        clear M;
                                    else
                                        %% Click on skel voxel or end point --> remove branch
                                        M = ((M == LUT(6))|(M == LUT(7)));
                                        M = regionGrow26(M,y,x,z,0,Inf);
                                        buf = imviewer.mask;
                                        buf(M>0) = LUT(2);
                                        setMask(imviewer,buf);
                                        clear M;
                                        clear buf;
                                    end
                                else
                                    disp('Seed voxel is null');
                                end
                            end
                            
                        case 'm' %% Toggle mask on/off
                            
                            ShowAnnot = ~ShowAnnot;
                            if ShowAnnot == false
                                setAlpha(imviewer,0);
                            else
                                setAlpha(imviewer,Alpha);
                            end
                            
                        case 'z' %% Toggle local Z projection on/off
                            
                            projmode = ~projmode;
                            if projmode == 1
                                RunProjViewer = RunProj;
                            else
                                RunProjViewer = 0;
                            end
                        
                        case 'v' %% 3D render
                
                            LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle);    
                            
                        case 'u' %% Undo
                            
                            if exist('old','var')
                                setMask(imviewer,uint8(old));
                                clear old;
                            end
                            
                        end
                 end
                    
            case 'f' %% Fit by gradient directed edge shrinking
                
               old = imviewer.mask;
               M = GradShrinkMask(old,I,2,5);
               imviewer.mask = M;
               disp('Done!');  
               
            case 'c' %% Classify 
                
                answer = questdlg('Classify? This will erase existing annotations');
                
                if answer(1) == 'Y'   
                
                %% Set state and reset annotations    
                Classify = true;
                ClassifierLoaded = 0;
                setMask(imviewer,uint8(zeros(size(imviewer.mask))));
                setAlpha(imviewer,Alpha);
                Lbl = imviewer.mask;
                
                %% Compute block features
                disp('Computing image block features...');
                h = waitbar(1,'Computing block features');
                
                %% Initialization
                nfeats_perscale = 2*(1+4+2*(Mode2D==0))*NShifts;
                if Mode2D == 0
                    shfs = [1 -1 1 -1 1 -1];
                    dims = [1 1 2 2 3 3];
                else
                    shfs = [1 -1 1 -1];
                    dims = [1 1 2 2];
                end
                Features = single(zeros(numel(I)/(BlckSizeXY*BlckSizeXY*BlckSizeZ),nfeats_perscale*numel(Scales)+2+(Mode2D==0)));
                
                %% Loop over scales
                for i = 1:numel(Scales)
                    
                    %% Downscaling
                    Is = I;
                    if Scales(i)>1
                        if Mode2D == 0
                            Is = resize3D(Is, [size(Is,1)/Scales(i) size(Is,2)/Scales(i) size(Is,3)/Scales(i)], 'bilinear');
                        else
                            Is = resize3D(Is, [size(Is,1)/Scales(i) size(Is,2)/Scales(i)], 'bilinear');
                        end
                    end
                    
                    %% Compute block mean
                    sumrows = sum(reshape(single(Is),BlckSizeXY,[]),1);
                    sumcols = sum(reshape(sumrows,size(Is,1)/BlckSizeXY,BlckSizeXY,[]),2);
                    clear sumrows;
                    I1 = reshape(sum(reshape(sumcols,size(Is,1)*size(Is,2)/(BlckSizeXY*BlckSizeXY),BlckSizeZ,[]),2),size(Is,1)/BlckSizeXY,size(Is,2)/BlckSizeXY,size(Is,3)/BlckSizeZ,1);
                    clear sumcols;
                    
                    %% Upscaling mean features
                    if Scales(i)>1
                        if Mode2D == 0
                            I1 = resize3D(I1, [size(I1,1)*Scales(i) size(I1,2)*Scales(i) size(I1,3)*Scales(i)], 'bilinear');
                        else
                            I1 = resize3D(I1, [size(I1,1)*Scales(i) size(I1,2)*Scales(i)], 'bilinear');
                        end
                    end
                    
                    %% Store mean features
                    Features(:,1+(i-1)*nfeats_perscale) = I1(:)/(BlckSizeXY*BlckSizeXY*BlckSizeZ);
                    for j = 1:4+2*(Mode2D==0)
                        I2 = circshift(I1,shfs(j),dims(j));
                        Features(:,1+j+(i-1)*nfeats_perscale) = I2(:)/(BlckSizeXY*BlckSizeXY*BlckSizeZ);
                    end
                    
                    %% Compute block std
                    sumrows = sum(reshape(single(Is).^2,BlckSizeXY,[]),1);
                    sumcols = sum(reshape(sumrows,size(Is,1)/BlckSizeXY,BlckSizeXY,[]),2);
                    clear sumrows;
                    I1 = reshape(sum(reshape(sumcols,size(Is,1)*size(Is,2)/(BlckSizeXY*BlckSizeXY),BlckSizeZ,[]),2),size(Is,1)/BlckSizeXY,size(Is,2)/BlckSizeXY,size(Is,3)/BlckSizeZ,1);
                    clear sumcols;
                    
                    %% Upscaling std features
                    if Scales(i)>1
                        if Mode2D == 0
                            I1 = resize3D(I1, [size(I1,1)*Scales(i) size(I1,2)*Scales(i) size(I1,3)*Scales(i)], 'bilinear');
                        else
                            I1 = resize3D(I1, [size(I1,1)*Scales(i) size(I1,2)*Scales(i)], 'bilinear');
                        end
                    end
                    
                    %% Store std features
                    Features(:,8+(i-1)*nfeats_perscale) = I1(:)/(BlckSizeXY*BlckSizeXY*BlckSizeZ);
                    Features(:,8+(i-1)*nfeats_perscale) = Features(:,2+(i-1)*nfeats_perscale)-Features(:,1+(i-1)*nfeats_perscale).^2;
                    for j = 1:4+2*(Mode2D==0)
                        I2 = circshift(I1,shfs(j),dims(j));
                        Features(:,8+j+(i-1)*nfeats_perscale) = I2(:)/(BlckSizeXY*BlckSizeXY*BlckSizeZ);
                        Features(:,8+j+(i-1)*nfeats_perscale) = Features(:,8+j+(i-1)*nfeats_perscale)-Features(:,1+(i-1)*nfeats_perscale).^2;
                    end
                end
                          
                %% Compute block maximum (for local classification)
                maxrows = max(reshape(I,BlckSizeXY,[]),[],1);
                maxcols = max(reshape(maxrows,size(I,1)/BlckSizeXY,BlckSizeXY,[]),[],2);clear maxrows;
                clear maxrows;
                I1 = reshape(max(reshape(maxcols,size(I,1)*size(I,2)/(BlckSizeXY*BlckSizeXY),BlckSizeZ,[]),[],2),size(I,1)/BlckSizeXY,size(I,2)/BlckSizeXY,size(I,3)/BlckSizeZ,1);
                clear maxcols;
                %% Upscaling block max for classification
                LocMax = resize3D(I1, size(I), 'bilinear');
                
                %% Compute spatial features
                if SpaFeat == 1
                    [X,Y,Z] = meshgrid(1:size(I1,2),1:size(I1,1),1:size(I1,3));
                    X = single(X);Y = single(Y);Z = single(Z);
                    clear I1;
                    Features = [Features Y(:) X(:) Z(:)];
                end
                
                %% Display blocks and features information
                disp([num2str(size(Features,1)) ' blocks']);
                disp([num2str(size(Features,2)) ' features']);
                if ishandle(h)
                	close(h);
                end
                
                %% Main loop (user interactions)
                fullpred = 1;
                firstiter = 1;
                while Classify == true
                    
                    %% Show color coded annotations
                    if ShowAnnot == true
                        setAlpha(imviewer,Alpha);
                    else
                        setAlpha(imviewer,0);
                    end
                    
                    %% Are there between 2 and 3 classes or is a classifier loaded?
                    LblsIn = Lbl(Lbl>0);
                    LblsIn = unique(LblsIn);
                    if (numel(LblsIn) >= 2 && numel(LblsIn) <= 3) || (ClassifierLoaded == 1)
                        
                        if (numel(LblsIn) >= 2 && numel(LblsIn) <= 3) && (ClassifierLoaded == 0)
                            
                            %% Compute block labels
                            maxrows = max(reshape(Lbl,BlckSizeXY,[]),[],1);
                            maxcols = max(reshape(maxrows,size(Lbl,1)/BlckSizeXY,BlckSizeXY,[]),[],2);clear maxrows;
                            clear maxrows;
                            I1 = reshape(max(reshape(maxcols,size(Lbl,1)*size(Lbl,2)/(BlckSizeXY*BlckSizeXY),BlckSizeZ,[]),[],2),size(Lbl,1)/BlckSizeXY,size(Lbl,2)/BlckSizeXY,size(Lbl,3)/BlckSizeZ,1);
                            clear maxcols;
                            NBlocksPerSlice = size(I1,1)*size(I1,2);
                            Labels = I1(:);
  
                            %% Train random forest
                            indstrain = find(Labels>0);
                            Labels = Labels(indstrain);
                            %% Set Class 2 label to 0
                            Labels(Labels == LUT(2)) = 0;
                            Feats = Features(indstrain,1:end-(2+(Mode2D==0))*(SpaFeat==0));
                            disp(['Training random forest from ' num2str(size(Features,2)+(2+(Mode2D==0))*(SpaFeat==1)) ' features']);
                            RFStruct = train_RF(double(Feats),Labels,'ntrees',NTrees,'method','c');
                            
                            %% SVM training (deprecated)
                            %disp('Training (SVM)...');
                            %%svmStruct = svmtrain(Feats,Labels,'kernel_function','polynomial','tolkkt',0.01,'kktviolationlevel',0.01);
                            %svmStruct = fitcsvm(Feats,Labels);
                            
                        end
                        
                        %% Predict blocks with random forest
                        O = uint8(zeros(size(I)));
                        Preds = uint8(zeros(size(Features,1),1));
                        mx = max(Features(:,1)); 
                        
                        %% Predict all blocks
                        if fullpred == 1
                            SigFeatsInds = find(Features(:,1) >= BckThr);
                            disp(['Random forest prediction from ' num2str(size(Features,2)+(2+(Mode2D==0))*(SpaFeat==1)) ' features and for ' num2str(numel(SigFeatsInds)) ' blocks']);
                            if length(SigFeatsInds) > 2*predchunks
                                h = waitbar(0,'Random Forest prediction');
                            end
                            for c = 1:predchunks:length(SigFeatsInds)
                                if length(SigFeatsInds) > 2*predchunks
                                    waitbar(c/length(SigFeatsInds));
                                end
                                Preds(SigFeatsInds(c:min(c+predchunks-1,length(SigFeatsInds)))) = uint8(eval_RF(double(Features(SigFeatsInds(c:min(c+predchunks-1,length(SigFeatsInds))),1:end-(2+(Mode2D==0))*(SpaFeat==0))),RFStruct,'oobe','y'));
                            end
                            clear SigFeatsInds;
                        else %% Predict only visible blocks
                            z = floor((getCurrentSlice(imviewer)-1)/BlckSizeZ);
                            FeaturesSlice = Features(1+z*NBlocksPerSlice:(z+1)*NBlocksPerSlice,:);
                            disp(['Random forest prediction from ' num2str(size(Features,2)+(2+(Mode2D==0))*(SpaFeat==1)) ' features and for ' num2str(size(FeaturesSlice,1)) ' blocks']);
                            Preds(1+z*NBlocksPerSlice:(z+1)*NBlocksPerSlice) = uint8(eval_RF(double(FeaturesSlice), RFStruct, 'oobe', 'y'));
                        end
                        
                        %% SVM prediction (deprecated)
                        %if fullpred == 1
                        %    SigFeatsInds = find(Features(:,1) >= BckThr);
                        %    disp(['Predicting (SVM)... ' num2str(num2str(numel(SigFeatsInds))) ' blocks']);
                        %    if length(SigFeatsInds) > 4*predchunks
                        %        h = waitbar(0,'SVM prediction');
                        %    end
                        %    for c = 1:predchunks:length(SigFeatsInds)
                        %        if length(SigFeatsInds) > 4*predchunks
                        %            waitbar(c/length(SigFeatsInds));
                        %        end
                        %        Preds(SigFeatsInds(c:min(c+predchunks-1,length(SigFeatsInds)))) = uint8(predict(svmStruct,Features(SigFeatsInds(c:min(c+predchunks-1,length(SigFeatsInds))),:)));
                        %    end
                        %    clear SigFeatsInds; 
                        %else
                        %    z = floor((getCurrentSlice(imviewer)-1)/BlckSizeZ);
                        %    FeaturesSlice = Features(1+z*NBlocksPerSlice:(z+1)*NBlocksPerSlice,:);
                        %    disp(['Predicting (SVM)... ' num2str(size(FeaturesSlice,1)) ' blocks']);
                        %    Preds(1+z*NBlocksPerSlice:(z+1)*NBlocksPerSlice) = uint8(predict(svmStruct,FeaturesSlice));
                        %end
                        
                        %% Upsample results
                        O = reshape(Preds,[size(I,1)/BlckSizeXY size(I,2)/BlckSizeXY size(I,3)/BlckSizeZ]);
                        O = round(resize3D(O,size(I),'nearest'));
                        
                        %% Local classification
                        disp('Local classification... ');
                        O = O.*uint8(I>=(LocMax*ThrLvl));
                        
                        %% Median filtering
                        if fullpred == 1
                            for i = 1:size(O,3)
                                O(:,:,i) = medfilt2(O(:,:,i),[3 3]);
                            end
                        else
                            z = getCurrentSlice(imviewer);
                            for i = max(z-BlckSizeZ,1):min(z+BlckSizeZ,size(O,3))
                                O(:,:,i) = medfilt2(O(:,:,i),[3 3]);
                            end
                        end               
                        
                        %% Clear waitbar
                        if ishandle(h)
                            close(h);
                        end
                        
                        %% Display mask
                        setMask(imviewer,O);
                        setAlpha(imviewer,Alpha);
                        MaskMode = false;
                        disp('Done!');
                        
                    else
                        
                        if firstiter==0
                            uiwait(msgbox('Pixels from at least 2 classes should be annotated!'));
                        end
                        
                    end
                    
                    %% User actions
                    userlabel = true;
                    set(handle,'Name','Classify: (0-3)   (Return)Pred   (Space)PredSlice  (A)nnotation/Pred   (S)ave   (L)oad          ||          (M)ask    (Z)proj   (V)olume   (O)ptions   e(X)it','NumberTitle','off');   
                    firstiter = 0; 
                    
                    while (userlabel == true)  
                       
                        %% Show color coded annotations
                        if ShowAnnot == true
                            setAlpha(imviewer,Alpha);
                        else
                            setAlpha(imviewer,0);
                        end
                        
                        %% Wait for user input
                        KeyRead = '';
                        while isempty(KeyRead)
                            pause(0.05);
                            currkey = KeyRead;
                        end
                        
                        switch currkey
                            
                            case 'x' %% Exit
                                answer = questdlg('Exit classification?');
                                if answer(1) == 'Y' 
                                    userlabel = false;
                                    Classify = false;
                                end
                                
                            case 'z' %% Toggle local Z projection
                                projmode = ~projmode;
                                if projmode == 1
                                    RunProjViewer = RunProj;
                                else
                                    RunProjViewer = 0;
                                end
                            
                            case 'v' %% 3D render
                
                                LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle);      
                                  
                            case 't' %% 
                                Features = Features(:,1:2);
                                
                            case 'return' %% Train + predict all blocks
                                
                                userlabel = false;
                                fullpred = 1; 
                                
                            case 'space' %% Train + predict visible blocks
                                
                                userlabel = false;
                                fullpred = 0;
                                
                            case 'a' %% Switch prediction mask / annotations visibility
                                
                                MaskMode = ~MaskMode;
                                if MaskMode == true
                                    setMask(imviewer,Lbl);
                                    setAlpha(imviewer,Alpha);
                                else
                                    if exist('O','var')
                                        setMask(imviewer,O);
                                        setAlpha(imviewer,Alpha);
                                    end
                                end
                                
                            case 'm' %% Toggle mask on/off
                                ShowAnnot = ~ShowAnnot;
                                if ShowAnnot == false
                                    setAlpha(imviewer,0);
                                else
                                    setAlpha(imviewer,Alpha);
                                end
                                
                            case 's' %% Save RF classifier
                                
                                if exist('RFStruct','var')
                                    [file,lastpath] = uiputfile('*.mat','Save RF classifier');
                                    fname = [lastpath file];
                                    save(fname,'RFStruct');
                                    disp('RF classifier saved');
                                end 
                                
                            case 'l' %% Load RF classifier
                                
                                [file,lastpath] = uigetfile('*.mat','Load RF classifier');
                                fname = [lastpath file];
                                load(fname,'RFStruct');
                                ClassifierLoaded = 1;
                                Lbl = Lbl*0;
                                disp('RF classifier loaded');
                                
                            case 'o' %% Options
                                
                                defaultValues = {num2str(NTrees),num2str(BckThr),num2str(ThrLvl),num2str(SpaFeat),num2str(RunProj),num2str(Alpha),};
                                titleBar = 'Set options';
                                userPrompt = {'Number of trees (5-100)','Background Threshold','Block Local Threshold (-1 - 1)','Use spatial features (0/1)','Z projection depth (3-127 slices)','Mask opacity (0.1-1)'};
                                caUserInput = inputdlg(userPrompt, titleBar, 1, defaultValues);
                                if ~isempty(caUserInput)
                                    NTrees = boundvar(str2double(caUserInput{1}),5,100,25);
                                    BckThr = boundvar(str2double(caUserInput{2}),-65535,65535,0);
                                    ThrLvl = boundvar(str2double(caUserInput{3}),-1,1,0);
                                    SpaFeat = boundvar(str2double(caUserInput{4}),0,1,0); 
                                    RunProj = boundvar(str2double(caUserInput{5}),3,127,0);
                                    Alpha = boundvar(str2double(caUserInput{6}),0.1,1,1);
                                end
                                if RunProjViewer > 0
                                    RunProjViewer = RunProj;
                                end
                                
                            otherwise %% Add annotation to Class 1-3
                            
                                if (ClassifierLoaded == 0)
                                    
                                c = str2double(currkey);
                                if ~isempty(c)
                                    if c>=0 && c<= 3
                                        %% Display command
                                        if c>0
                                            if c == 1
                                                set(handle,'Name','Class background','NumberTitle','off');
                                            else
                                                set(handle,'Name',['Class' num2str(c)],'NumberTitle','off');
                                            end
                                        else
                                            set(handle,'Name','Erase label','NumberTitle','off'); 
                                        end
                                        
                                        %% User rectangular annotation
                                        h = drawrectangle();
                                        coords = round(h.Position);
                                        
                                        delete(h);
                                        coords(1) = max(1,coords(1));
                                        coords(2) = max(1,coords(2));
                                        z = round(getCurrentSlice(imviewer));

                                        %% Update labels
                                        slice = uint8(zeros(size(I,1),size(I,2)));
                                        old = Lbl(:,:,z);
                                        if c>0
                                            slice(coords(:,2):min(coords(:,2)+coords(:,4)-1,size(slice,1)),coords(:,1):min(coords(:,1)+coords(:,3)-1,size(slice,2))) = LUT(1+c);
                                            old(slice>0) = 0;
                                            Lbl(:,:,z) = old+slice;
                                        else
                                            slice(coords(:,2):min(coords(:,2)+coords(:,4)-1,size(slice,1)),coords(:,1):min(coords(:,1)+coords(:,3)-1,size(slice,2))) = 1;
                                            old(slice>0) = 0;
                                            Lbl(:,:,z) = old;
                                        end
                                        if MaskMode == true
                                            setMask(imviewer,Lbl);
                                        end

                                    end
                                end
                                
                                else
                                    
                                    uiwait(msgbox('Loaded classifier cannot be re-trained from annotations'));
                                
                                end
                        end
                        end
                end
            end
            if exist('O','var')
                setMask(imviewer,O);
            end
            clear Features;
            clear Lbl;
            clear O;
            clear I1;
            clear Preds;
            
        end
        
        %% Exception handling (prevent crash)
        catch ME
           if strcmp(ME.message,'Stopped by user')
               exit;
           else
               warning('Something went wrong, operation canceled');
               disp(ME.message);
           end
        end
        
    end
       
end
close all;

function keyPress(src, e)
    global KeyRead;
    KeyRead = e.Key;
end