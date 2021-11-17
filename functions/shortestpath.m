function inds = shortestpath(Dmap,y,x,z,y2,x2,z2)

    Dmap(1,:,:) = Inf;Dmap(end,:,:) = Inf;Dmap(:,1,:) = Inf;Dmap(:,end,:) = Inf;Dmap(:,:,1) = Inf;Dmap(:,:,end) = Inf;
    %Pre-allocate buffer
    inds = zeros(65535,1);
    inds(1) = sub2ind(size(Dmap),y,x,z);
    curind = inds(1);
    finalind = sub2ind(size(Dmap),y2,x2,z2);
    
    %% Neighbors indexing
    ncol = size(Dmap,1);
    nslice = size(Dmap,1)*size(Dmap,2);
    n = [1 -1 ncol -ncol 1+ncol 1-ncol -1+ncol -1-ncol...
         nslice nslice+1 nslice-1 nslice+ncol nslice-ncol nslice+1+ncol nslice+1-ncol nslice-1+ncol nslice-1-ncol...
         -nslice -nslice+1 -nslice-1 -nslice+ncol -nslice-ncol -nslice+1+ncol -nslice+1-ncol -nslice-1+ncol -nslice-1-ncol];
    N = size(Dmap,1)*size(Dmap,2)*size(Dmap,3);
    
    %% Iterate until completion
    it = 1;
    while curind ~= finalind && it < 65535
        Cand = n + curind;
        Cand(Cand<1) = [];
        Cand(Cand>N) = [];
        [mn, mnind] = min(Dmap(Cand));
        if mn < Inf
            curind = Cand(mnind);
            it = it + 1;
            inds(it) = curind;
        else
            it = 65535;
            inds = [];
            disp('No connecting paths');
        end
    end
    
    if ~isempty(inds)
        inds = inds(1:it);
    end
    
end

