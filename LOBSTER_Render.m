function LOBSTER_Render(imviewer, I, LUT, LUT3d, NCols, PointSize, ZSpan, ZRatio, DownSamp, handle)

    %% 6 connectivity structure element
    se = uint8(zeros(3,3,3));
    se(:,:,1) = [0 0 0;0 1 0;0 0 0];
    se(:,:,2) = [0 1 0;1 1 1;0 1 0];
    se(:,:,3) = [0 0 0;0 1 0;0 0 0];

    %% Set rendering mode depending on mask classes
    mx = max(imviewer.mask(:));
    RenderPointSize = PointSize;
    if mx == LUT(2)
        disp('3D rendering: Class 1 objects)');
        ClearClass1 = 0;
        ConnAnalysis = 1; 
    elseif mx == LUT(3)
        disp('3D rendering:Class 2 objects (Class 1 transparent)');
        ClearClass1 = 1;
        ConnAnalysis = 2;
    elseif mx == LUT(4)
        disp('3D rendering: Class 2 and 3 (Class 1 transparent)');
        ClearClass1 = 1;
        ConnAnalysis = 0;
    else
        disp('3D rendering: Skeleton mode (Class 1 transparent)');
        ClearClass1 = 1;
        ConnAnalysis = 0;
    end

    %% Get user input for render region
    figure(handle);
    h = drawrectangle();
    pos = round(h.Position);
    if (pos(3)==0)&&(pos(4)==0)
        pos(1) = pos(1)-32;
        pos(2) = pos(2)-32;
        pos(3) = max(pos(3),64);
        pos(4) = max(pos(3),64);
    end
        
    %% Convert image to 8 bit if necessary
    if isa(I, 'uint16')
        scl = 255/single(max(I(:)));
        I = uint8(single(I)*scl);
    end

    %% Crop and convert to 8-bit for rendering
    if isinf(ZSpan)
        Icrop = I(max(floor(pos(2)),1):min(floor(pos(2)+pos(4)),size(I,1)),max(floor(pos(1)),1):min(floor(pos(1)+pos(3)),size(I,2)),1:end);
        Ocrop = imviewer.mask(max(floor(pos(2)),1):min(floor(pos(2)+pos(4)),size(I,1)),max(floor(pos(1)),1):min(floor(pos(1)+pos(3)),size(I,2)),1:end);
    else
        Icrop = I(max(floor(pos(2)),1):min(floor(pos(2)+pos(4)),size(I,1)),max(floor(pos(1)),1):min(floor(pos(1)+pos(3)),size(I,2)),max(getCurrentSlice(imviewer)-floor(ZSpan/2)+1,1):min(getCurrentSlice(imviewer)+floor(ZSpan/2),size(I,3)));
        Ocrop = imviewer.mask(max(floor(pos(2)),1):min(floor(pos(2)+pos(4)),size(I,1)),max(floor(pos(1)),1):min(floor(pos(1)+pos(3)),size(I,2)),max(getCurrentSlice(imviewer)-floor(ZSpan/2)+1,1):min(getCurrentSlice(imviewer)+floor(ZSpan/2),size(I,3)));
    end

    %% Compute required XY padding to pad to cube
    maxXYZ = max([ceil(size(Icrop,1)/8)*8 ceil(size(Icrop,2)/8)*8 ceil(size(Icrop,3)/8)*8]);
    CubeRatio = maxXYZ/size(Icrop,3);

    %% Pad to cube for rendering
    Icrop = permute(padarray(Icrop,[maxXYZ-size(Icrop,1) maxXYZ-size(Icrop,2) max(8-size(Icrop,3),0)],0,'pre'),[2 1 3]);
    Ocrop = permute(padarray(Ocrop,[maxXYZ-size(Ocrop,1) maxXYZ-size(Ocrop,2) max(8-maxXYZ-size(Ocrop,3),0)],0,'pre'),[2 1 3]);

    %% Compute effective ZRatio (ZScale)
    if(exist('ZRatio','var'))
        ZScale = (ZRatio/DownSamp)/CubeRatio;
        disp(['ZRatio: ' num2str(ZRatio)]);
        disp(['Rendered volume: ' num2str(size(Icrop,2)) ' x ' num2str(size(Icrop,1)) ' x ' num2str(size(Icrop,3))]);
    end
    disp(['XY Down Sampling: ' num2str(DownSamp)]);
    disp(['Rendering ZScale: ' num2str(ZScale)]);

    %% Make class 1 transparent
    if ClearClass1
        Ocrop(Ocrop==0) = LUT(2);
    end

    %% Find edges
    Ocropshrink = imerode(Ocrop,se);
    Edges = Ocrop-Ocropshrink;

    %% Compute all surface points and assign color
    if numel(Edges)>0
        ObjIndx = find(Edges>0);	
        [i,j,k] = ind2sub(size(Ocrop),ObjIndx);
        FV.vertices = [i j CubeRatio*k];
        FV.vertices = FV.vertices*(2/max(size(Icrop)))-1;
        FV.vertices(:,3) = FV.vertices(:,3)*ZScale;
        FV.faces = repmat((1:numel(i)).',1,3);
        if ConnAnalysis >= 1
            Lcrop = bwlabeln(Ocrop==LUT(1+ConnAnalysis));
            cols = uint8(LUT3d(1+mod(Lcrop(ObjIndx),NCols),:));
        else
            cols = uint8(LUT3d(1+mod(Ocrop(ObjIndx),NCols),:));
        end
        FV.colors = cols;
        FV.pointsize = RenderPointSize;
    end
    FV.points = 1;
    FV.ZScale = ZScale;

    %% Mesh mode  
    %FV = isosurface(Ocrop,LUT(3)-1);
    %FV.Normals = isonormals(Ocrop,FV.vertices);
    %FV.vertices(:,3) = FV.vertices(:,3)*CubeRatio;
    %FV.vertices = FV.vertices*(2/max(size(Icrop)))-1;
    %FV.vertices(:,3) = FV.vertices(:,3)*ZScale;
    %buf = FV.vertices(:,1);
    %FV.vertices(:,1) = FV.vertices(:,2);
    %FV.vertices(:,2) = buf; 
    %FV.ZScale = ZScale;
    %FV.pointsize = PointSize;
    %LUT3d = round(linspecer(NCols)*255);
    %rng(0);
    %LUT3d = LUT3d(randperm(NCols),:);
    %FV.colors = uint8(repmat(LUT3d(1,:),numel(FV.vertices),1));

    %% 3D rendering
    LoadTaoOpenGl;
    Form1 = Init_OpenGL_Window('new','render3d_and_mesh',Icrop, FV);
    Form1.WindowState = System.Windows.Forms.FormWindowState.Normal;
    Form1.TopMost = true;
    while Form1.Visible == 1
        pause(0.05); 
    end
    delete(h);