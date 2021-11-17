function [O] = fxg_sLoG3DLocMax3D(I, params)

    % Apply 3D LoG + invert + detect 3D local intensity maxima (mark seeds).
    %
    % Input: 3D grayscale image
    % Output: 3D seed mask
    %
    % Sample journal: <a href="matlab:JENI('CellPilar3D_DetLog3DLocMax3D.jls');">CellPilar3D_DetLog3DLocMax3D.jls</a>
    %
    % Parameters:
    % Sigmas:           X,Y,Z sigmas of Gaussian blur pre-filter (vector, pix)
    % LocalMaxBox:      Local intensity maxima search box size (pix)     
    % MinLoG:           Minimum LoG response for object detection (normalized to 1)

    %% Parameters
    Sigmas = params.Sigmas;
    LocalMaxBox = params.LocalMaxBox;
    MinLoG = params.MinLoG;

    if ~isempty(I)

        %% Initialize filters
        LocMxse = ones(LocalMaxBox(1),LocalMaxBox(2),LocalMaxBox(3));
        LocMxse(ceil(end/2),ceil(end/2),ceil(end/2)) = 0; 
        
        disp('Filtering stack...');
        %% 3DLoG filter
        I = single(I);
        G2 = fspecial('gauss',[round(5*Sigmas(2)) 1], Sigmas(2));
        G1 = fspecial('gauss',[round(5*Sigmas(1)) 1], Sigmas(1));
        G3 = fspecial('gauss',[round(5*Sigmas(3)) 1], Sigmas(3));
        If = imfilter(I,G2,'same','symmetric'); 
        If = imfilter(If,permute(G1,[2 1 3]),'same','symmetric'); 
        If = imfilter(If,permute(G3,[3 2 1]),'same','symmetric'); 
        [Gx,Gy,Gz] = gradient(If);clear If;
        [Gxx,Gxy,Gxz] = gradient(Gx);clear Gxy;clear Gxz;clear Gx;
        [Gyx,Gyy,Gyz] = gradient(Gy);clear Gy;clear Gyx;clear Gyz;Gxx = Gxx+Gyy;clear Gyy;
        [Gzx,Gzy,Gzz] = gradient(Gz);clear Gz;clear Gzx;clear Gzy;Gxx = Gxx+Gzz;clear Gzz;
        ILog = -Gxx;clear Gxx;clear Gyy;clear Gzz;

        %% Normalize LoG response
        ILog = ILog/max(ILog(:));
        
        %% LoG local maxima detection
        disp('Detecting local maxima...');
        ILogMax = imdilate(ILog,LocMxse);
        O = (ILog > ILogMax)&(ILog >= MinLoG);
        
        %% Create seed mask
        Seeds = fxm_sMarkObjCentroids(O,[]);
        %O = uint8(100*imdilate(O,ones(5,5)));
        
        O = uint8(zeros(size(O)));
        O(find(Seeds)) = 200;
        
    else
        
    	O = [];
        
    end
    
end