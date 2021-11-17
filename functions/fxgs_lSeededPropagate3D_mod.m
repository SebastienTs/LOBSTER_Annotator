function [L] = fxgs_lSeededPropagate3D_mod(I,L,M,params)
    
    % Propagate seeds through low intensity variation regions.
    %
    % Sample journal: <a href="matlab:JENI('CellPilar3D_LogLocMaxLocThrPropagate3D.jls');">CellPilar3D_LogLocMaxLocThrPropagate3D.jls</a>
    %
    % Input: 3D original image, 3D binary or label mask, optional 3D binary mask
    % Output: 3D label mask
    %
    % Parameters:
    % Power:                Apply power law to intensity image prior to processing
    % AnalyzeCC:            Set to 0 for input label mask is passed, 1 for binary mask
    
    %% Parameters
    Power = params.Power;
    AnalyzeCC = params.AnalyzeCC;

    if ~isempty(I)
        
		if AnalyzeCC
            L = single(bwlabeln(L,6));
        else
            L = single(L);
        end
        
        %% Set mask borders to 0 since no image border check is performed in Propgate_3D
        M(1,:,:) = 0;M(end,:,:) = 0;M(:,1,:) = 0;M(:,end,:) = 0;M(:,:,1) = 0;M(:,:,end) = 0;
        
        if Power ~= 1
            L = Propagate_3D_single(L,(I).^Power,logical(M),single(size(I,3)));
        else
            L = Propagate_3D_single(L,I,logical(M),single(size(I,3)));
        end
        L = reshape(L,size(I));
        
    else
        
        L = [];
    
    end
    
end