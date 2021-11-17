function M = LayTube(M,Points,Rad,fillval,ZRatio)

    %% Read mask size
    Dim = size(M);

    %% Compute ZRatio adjusted sphere voxel offsets
    [y,x,z] = meshgrid(-Rad:Rad,-Rad:Rad,-round(Rad/ZRatio):round(Rad/ZRatio));
    m = sqrt(x.^2 + y.^2 + (z*ZRatio).^2); 
    b = (m <= Rad); 
    [yp,xp,zp] = ind2sub(size(b),find(b==1));
    offs = yp-Rad-1+(xp-Rad-1)*Dim(1)+(zp-round(Rad/ZRatio)-1)*Dim(1)*Dim(2);

    %% Read point coordinates
    X = Points(:,2).';
    Y = Points(:,1).';
    Z = Points(:,3).';

    %% Compute spline distance and steps to densely lay spheres (4 per unit radius)
    dist = sqrt(sum(diff(X).^2+diff(Y).^2+diff((Z*ZRatio)).^2));
    Nsteps = ceil(4*dist/Rad);
    
    %% Compute spline function
    xyz = [X; Y; Z];
    f = cscvn(xyz(:,[1:end]));

    %% Plot control points and spline function
    %plot3(X,Y,Z,'x');
    %hold on;
    %fnplt(f,'r',2);
    
    %% Write spline cylinder to image mask
    t = f.breaks(end)/Nsteps*[0:Nsteps];
    coords = round(fnval(f,t));
    inds = coords(1,:)+coords(2,:)*size(M,1)+coords(3,:)*size(M,1)*size(M,2);
    fillinds = round(repmat(inds,numel(offs),1)+repmat(offs,1,numel(inds)));
    fillinds(fillinds<1) = [];fillinds(fillinds>numel(M)) = [];
    M(fillinds) = fillval;
    %figure;imagesc(max(M,[],3));
    disp('Done');
    
end