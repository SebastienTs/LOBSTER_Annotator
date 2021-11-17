function M = GradShrinkMask(M,I,GRad,MaxIter)

    %% Variables
    Sx = size(I,1);
    Sz = size(I,1)*size(I,2);
    nhood26 = ones(3,3,3);
    nhood6= zeros(3,3,3);
    nhood6(14) = 1;nhood6(13) = 1;nhood6(15) = 1;
    nhood6(11) = 1;nhood6(17) = 1;nhood6(5) = 1;nhood6(23) = 1;
    
    %% Filter image
    If = single(imgaussfilt3(I,GRad));
    
    %% Estimate intensity gradient
    [Ifx,Ify,Ifz] = gradient(If);
    GradMag = sqrt(Ifx.^2+Ify.^2+Ifz.^2);
    
    %% Cleanup, close small holes
    %M = imclose(M,nhood26);
    
    %% Remove mask border pixels
    M(1,:,:) = 0;M(end,:,:) = 0;
    M(:,1,:) = 0;M(:,end,:) = 0;
    M(:,:,1) = 0;M(:,:,end) = 0;
    
    %% Main loop
    for i = 1:MaxIter
        %% Detect Mask edges
        Edg = logical(M)-logical(imerode(M,nhood6));
        cands = find(Edg & ((~isnan(GradMag))&(~(GradMag==0))) );
		offs = sign(Ify(cands)) + Sx*sign(Ifx(cands)) + Sz*sign(Ifz(cands));
        vals = GradMag(cands)-GradMag(cands+offs);
        M(cands(vals<0)) = 0;
        %M = imopen(M,nhood6);
    end 
    
    %% Cleanup, close small holes
    %M = imclose(M,nhood26);
    
end