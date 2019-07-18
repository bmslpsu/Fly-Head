% make rotation pattern with 4 different widths, first 4 indices of Y
% contain the phi version and 5-8 are the rev_phi version

x_num = 96;  %% x contains rotations all around the display
y_num = 6;   %% y is the indices into the spatial freq. and phi\rev phi

Pats = zeros(1, 96, x_num, y_num);

% first the 4 phi pattern
Pats(:,:,1,1) = repmat([3*ones(1,12) , 6*ones(1,12)], 1, 4);
Pats(:,:,1,2) = repmat([3*ones(1,8) , 6*ones(1,8)], 1, 6);
Pats(:,:,1,3) = repmat([3*ones(1,4) 6*ones(1,4)], 1, 12);

% the 4 rev phi patterns
Pats(:,:,1,4) = repmat([3*ones(1,12) , 6*ones(1,12)], 1, 4);
Pats(:,:,1,5) = repmat([3*ones(1,8) , 6*ones(1,8)], 1, 6);
Pats(:,:,1,6) = repmat([3*ones(1,4) 6*ones(1,4)], 1, 12);

for j = 2:x_num 			%use ShiftMatrixPats to rotate stripe image
    for k = 1:y_num
        Pats(:,:,j,k) = ShiftMatrix(Pats(:,:,j-1,k),1,'r','y');
    end
end

% for the 4 rev phi pattern, on alternating frames replace the bright
% levels (positive 1) with the dim level (here -1)

for k = 4:6
    for j = 2:2:96 % accounting for each frame in x except frame 1
        Ind_1 = find(Pats(:,:,j,k) == 6);
        [I,J] = ind2sub([1 96],Ind_1);
        for m = 1:length(I)
            Pats(I(m),J(m),j,k) = 0;
        end
    end
    
end

% make a space-time plot of the stimuli
cmap = [0 0 0; 0 1/6 0;0 2/6 0; 0 3/6 0;0 4/6 0; 0 5/6 0;0 6/6 0];
colormap(cmap);
for k = 1:6
    figure(7)
    subplot(2,3,k)
    imagesc(squeeze(Pats(1,:,:,k))');colormap(cmap);
    hold on;
    ylabel('time')
    xlabel('space');
    if k == 1
        title('standard motion');
    elseif k == 4
     title('reverse-phi motion');
    end
end