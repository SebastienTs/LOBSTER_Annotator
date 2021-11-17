function [O] = fxg_sLoGLocMax(I, params)

    % Apply LoG + invert + detect 3D local intensity maxima (mark seeds).
    %
    % Sample journal: <a href="matlab:JENI('NucleiNoisy_DetLogLocMax.jl');">NucleiNoisy_DetLogLocMax.jl</a>
    %
    % Input: 2D grayscale image
    % Output: 2D seed mask
    %
    % Parameters:
    % GRad:             Gaussian blur pre-filter radius (pix)
    % LocalMaxBox:      Size of local intensity maxima search box (pix)     
    % MinLoG:           Minimum LoG response for object detection (normalized to 1)

    % Parameters
    GRad = params.GRad;
    LocalMaxBox = params.LocalMaxBox;
    MinLoG = params.MinLoG;
    
    if ~isempty(I)

        %% Initialize filters
        H = fspecial('log',3*GRad,GRad);
        LocMxse = ones(LocalMaxBox(1),LocalMaxBox(2));
        LocMxse(ceil(end/2),ceil(end/2)) = 0;

        %% XY LoG filter
        I = double(I);
        ILog = -imfilter(I,H,'same','symmetric');

        %% Normalize LoG response
        ILog = ILog/max(ILog(:));
        
        %% Local maxima detection
        Maxima = ILog > imdilate(ILog,LocMxse);
        Maxima = Maxima.*(ILog>MinLoG);
        
        %% Local maxima filtering (intensity)
        Seeds = fxm_sMarkObjCentroids(Maxima,[]);
        O = uint8(100*Maxima);
        O(find(Seeds)) = 200;

    else
        
        O = [];
        
    end
    
end