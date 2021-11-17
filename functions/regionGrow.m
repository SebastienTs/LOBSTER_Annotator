function M = regionGrow(I,y,x,z,tol,mode,dst)

    %% Initial mask
    M = uint8(zeros(size(I)));

    if y > 1 && x > 1 && y < size(I,1) && x < size(I,2)

        %% Neighbors indexing
        n = [1 -1 size(I,1) -size(I,1) size(I,1)*size(I,2) -size(I,1)*size(I,2)];
        N = size(I,1)*size(I,2)*size(I,3);
        
        %% Growth mask: prevent x/y edges to be conquered, initialize seed
        M(1,:,:) = 1;M(end,:,:) = 1;M(:,1,:) = 1;M(:,end,:) = 1;
        M(y,x,z) = 1;
        Cur = sub2ind(size(M),y,x,z);
        Thrmin = I(y,x,z)*(1-tol);
        Thrmax = I(y,x,z)*(1+tol);
        
        %% Iterate until completion
        it = 0;
        while ~isempty(Cur) && it<dst
            Cand = repmat(Cur,1,numel(n))+repmat(n,numel(Cur),1);
            Cand = Cand(:);
            Cand = Cand((Cand>0)&(Cand<N));
            if mode == 1
                Cand = Cand((M(Cand) == 0)&((I(Cand) >= Thrmin)));
            else
                Cand = Cand((M(Cand) == 0)&((I(Cand) >= Thrmin))&((I(Cand) <= Thrmax)));
            end
            Cand = unique(Cand);
            M(Cand) = 1;
            Cand = Cand();
            Cur = Cand;
            it = it + 1;
        end

        %% Remove x/y edge border
        M(1,:,:) = 0;M(end,:,:) = 0;M(:,1,:) = 0;M(:,end,:) = 0;
        
    end
    
end