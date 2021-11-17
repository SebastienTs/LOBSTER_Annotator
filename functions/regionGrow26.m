function M = regionGrow26(I,y,x,z,tol,dst)

    M = uint8(zeros(size(I)));

    if y > 1 && x > 1 && y < size(I,1) && x < size(I,2)

        %% Neighbors indexing
        a = size(I,1);
        b = size(I,2);
        n = [1 -1 a -a -1+a -1-a 1+a 1-a -a*b -a*b+1 -a*b-1 -a*b+a -a*b-a -a*b-1+a -a*b-1-a -a*b+1+a -a*b+1-a -a*b a*b a*b+1 a*b-1 a*b+a a*b-a a*b-1+a a*b-1-a a*b+1+a a*b+1-a];
        N = size(I,1)*size(I,2)*size(I,3);
        
        %% Growth mask: prevent x/y edges to be conquered, initialize seed
        M(1,:,:) = 1;M(end,:,:) = 1;M(:,1,:) = 1;M(:,end,:) = 1;
        M(y,x,z) = 1;
        Cur = sub2ind(size(M),y,x,z);
        
        if I(y,x,z) > 0
            Thrmin = I(y,x,z)*(1-tol);
            Thrmax = I(y,x,z)*(1+tol);

            %% Iterate until completion
            it = 0;
            while ~isempty(Cur) && it<dst
                Cand = repmat(Cur,1,numel(n))+repmat(n,numel(Cur),1);
                Cand = Cand(:);
                Cand = Cand((Cand>0)&(Cand<N));
                %Cand = Cand((M(Cand) == 0)&((I(Cand) >= Thrmin)));
                Cand = Cand((M(Cand) == 0)&((I(Cand) >= Thrmin))&((I(Cand) <= Thrmax)));
                Cand = unique(Cand);
                M(Cand) = 1;
                Cand = Cand();
                Cur = Cand;
                it = it + 1;
            end
        else
            disp('Seed voxel is null');
        end

        %% Remove x/y edge border
        M(1,:,:) = 0;M(end,:,:) = 0;M(:,1,:) = 0;M(:,end,:) = 0;
        
    end
    
end


