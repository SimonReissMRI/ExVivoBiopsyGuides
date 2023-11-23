%%

addpath('matlab_functions')

%% Load Multi-Echo DICOM Dataset
[filename,pathname]=uigetfile('*.dcm','Select all magnitude images','MultiSelect','On');

no_images = size(filename,2);

no_echoes = 3;

imsize = size(dicomread(fullfile(pathname,char(filename(1))))).*1;
imsize(3) = size(filename,2)/3;
images = zeros(imsize(1),imsize(2),imsize(3),no_echoes);
t = zeros(no_echoes,1);

temp = dicominfo(fullfile(pathname,char(filename(1))));
t(1) = temp.EchoTime;
temp = dicominfo(fullfile(pathname,char(filename(1+no_images/3))));
t(2) = temp.EchoTime;
temp = dicominfo(fullfile(pathname,char(filename(1+2*no_images/3))));
t(3) = temp.EchoTime;

dx = temp.PixelSpacing(1);
dy = temp.PixelSpacing(2);
dz = temp.SliceThickness;

for k=1:imsize(3)
    
    images(:,:,k,1)  = imresize(mean(dicomread(fullfile(pathname,char(filename(k)))),3),1);   
    images(:,:,k,2)  = imresize(mean(dicomread(fullfile(pathname,char(filename(k+no_images/3)))),3),1);  
    images(:,:,k,3)  = imresize(mean(dicomread(fullfile(pathname,char(filename(k+2*no_images/3)))),3),1);  
    
    clear temp

end

figure(1)
imagesc(images(:,:,end/2,1))
daspect([1 1 1])
colormap gray

images_TE1 = squeeze(images(:,:,:,1));

clear filename
%% Calculate R2* Map

T2starmap = zeros(imsize(1),imsize(2),imsize(3));
X = [ones(length(t),1) t];

log_images = log(images);

for i = 1:imsize(1)
   for j = 1:imsize(2)
       parfor k = 1:imsize(3)
            data = double(squeeze(log_images(i,j,k,:)));
            b = X\data;           
            T2starmap(i,j,k) = -1/b(2); 
       end
     
   end
end

R2starmap = 1000./abs(T2starmap);

clear T2starmap images log_images 

%% Import R2* Map

% load('R2star_map.mat')

%% Automatic Thresholding

% bin the data according to signal values
[N, edges] = histcounts(images_TE1(:));
% get center of signal bins
signalVals = diff(edges)/2 + edges(1:end-1);
% obtain some starting parameters
[pks, locs] = findpeaks(N, 'MinPeakDistance', 0.8*numel(N)/2);
% sigma of rician (position of first peak)
b = signalVals(locs(1));
% amplitude-factor of rician (maximum of rician)
a = pks(1) * sqrt(2)*exp(1/2)*b;
% for noise, the mean value nu is 0 and the rician simplifies to:
fitfcn = fittype('a*x/b^2*exp(-(x^2)/(2*b^2))');
f1 = fit(signalVals', N', fitfcn, 'start', [a b], ...
                                'lower', [0 signalVals(locs(1))/4], ...
                                'upper', [Inf max(signalVals)]);

mag_thresh = 5*f1.b;

%% Set Manual Magnitude Threshold

% figure(2)
% histogram(images_TE1(:), 180,'FaceColor', [0.63 0.73 0.93], 'EdgeColor', 'k');
% set(gca, 'YScale', 'log');
% xlabel('Signal Intensity [a.u.]');
% ylabel('Number of Voxels');
% xlim([0 1300]);
% xticks(0:200:1300);
% ylim([1.9e4 1e8]);
% 
% prompt = 'What is the treshold? ';
% mag_thresh = input(prompt);

%% Magnitude Tresholding
mag_mask_raw = ones(size(images_TE1));
mag_mask_raw(images_TE1<mag_thresh)=0;

mag_mask_inter = imresize3(imresize3(mag_mask_raw,1/4),4);
mag_mask = RegGrow(mag_mask_inter,0.1,[1/2*size(mag_mask_inter,1) 1/2*size(mag_mask_inter,2) 1/2*size(mag_mask_inter,3)]);

figure(3)
subplot(1,3,1)
ax = gca;
imagesc(squeeze(mag_mask(:,:,end/2)));
set(gca,'box','off');
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';
daspect([1 1 1]);

subplot(1,3,2)
ax = gca;
imagesc(squeeze(mag_mask(:,end/2,:)));
daspect([1 2 1]);
set(gca,'box','off');
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';

subplot(1,3,3)
ax = gca;
imagesc(squeeze(mag_mask(end/2,:,:)));
daspect([1 2 1]);
set(gca,'box','off');
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';


figure(4)
subplot(1,3,1)
imagesc(squeeze(mag_mask(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mag_mask(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mag_mask(end/2,:,:)));
daspect([1 2 1]);

figure(5)
subplot(1,3,1)
imagesc(squeeze(mag_mask(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mag_mask(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mag_mask(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);


%% R2* Tresholding

% bin the data according to signal values
[N, edges] = histcounts(R2starmap(:));
% get center of signal bins
signalVals = diff(edges)/2 + edges(1:end-1);
% obtain some starting parameters
[pks, locs] = findpeaks(N, 'MinPeakDistance', 0.8*numel(N)/2);
% sigma of rician (position of first peak)
b = signalVals(locs(1));
% amplitude-factor of rician (maximum of rician)
a = pks(1) * sqrt(2)*exp(1/2)*b;
% for noise, the mean value nu is 0 and the rician simplifies to:
fitfcn = fittype('a*x/b^2*exp(-(x^2)/(2*b^2))');
f1 = fit(signalVals', N', fitfcn, 'start', [a b], ...
                                'lower', [0 signalVals(locs(1))/4], ...
                                'upper', [Inf max(signalVals)]);

R2star_thresh = 3*f1.b;

R2star_mask = ones(size(R2starmap));
R2star_mask(R2starmap<R2star_thresh)=0;


figure(6)
subplot(1,3,1)
imagesc(squeeze(R2star_mask(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(R2star_mask(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(R2star_mask(end/2,:,:)));
daspect([1 2 1]);

figure(7)
subplot(1,3,1)
imagesc(squeeze(R2star_mask(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(R2star_mask(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(R2star_mask(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);

%% Create Combined Mask

mask_comb = R2star_mask.*mag_mask;


figure(8)
subplot(1,3,1)
imagesc(squeeze(mask_comb(:,:,end/2)).*squeeze(mag_mask(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mask_comb(:,end/2,:)).*squeeze(mag_mask(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mask_comb(end/2,:,:)).*squeeze(mag_mask(end/2,:,:)));
daspect([1 2 1]);

figure(8)
subplot(1,3,1)
imagesc(squeeze(mask_comb(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mask_comb(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mask_comb(end/2,:,:)));
daspect([1 2 1]);

figure(9)
subplot(1,3,1)
imagesc(squeeze(mask_comb(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mask_comb(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mask_comb(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);



%%  Remove Small Loose Parts

fprintf('Removing loose parts...\n');
cc = bwconncomp(mask_comb);
% find largest connected component
[~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
%rebuild component
mask_final = zeros(size(mask_comb));
mask_final(cc.PixelIdxList{index}) = 1;

for i = 1:size(mask_final,3)
    cc = bwconncomp(squeeze(mask_final(:,:,i)));
    % find largest connected component
    [~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
    %rebuild component
    temp = zeros(size(mask_final,1),size(mask_final,2));
    if (~isempty(index))
        temp(cc.PixelIdxList{index}) = 1;
    end
    mask_final(:,:,i) = temp.*squeeze(mask_final(:,:,i));
end

for i = 1:size(mask_final,1)
    cc = bwconncomp(squeeze(mask_final(i,:,:)));
    % find largest connected component
    [~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
    %rebuild component
    temp = zeros(size(mask_final,2),size(mask_final,3));
    if (~isempty(index))
        temp(cc.PixelIdxList{index}) = 1;
    end
    mask_final(i,:,:) = temp.*squeeze(mask_final(i,:,:));
end

for i = 1:size(mask_final,2)
    cc = bwconncomp(squeeze(mask_final(:,i,:)));
    % find largest connected component
    [~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
    %rebuild component
    temp = zeros(size(mask_final,1),size(mask_final,3));
    if (~isempty(index))
        temp(cc.PixelIdxList{index}) = 1;
    end
    mask_final(:,i,:) = temp.*squeeze(mask_final(:,i,:));
end

figure(10)
subplot(1,3,1)
imagesc(squeeze(mask_final(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mask_final(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mask_final(end/2,:,:)));
daspect([1 2 1]);

figure(11)
subplot(1,3,1)
imagesc(squeeze(mask_final(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(mask_final(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(mask_final(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);




%% Cut top and Connect

cut_slice = 35;

if (mod(cut_slice,2)==0)
    cut_slice = cut_slice+1;
end    

mask_final_cut = mask_final;
mask_final_cut(:,:,1:cut_slice)=0; 

fprintf('connect...\n');
cc2 = bwconncomp(1-mask_final_cut,6);
% find largest connected component
[~,index] = max(cellfun(@(x) numel(x),cc2.PixelIdxList));
%rebuild component
inv = zeros(size(mask_final_cut));
inv(cc2.PixelIdxList{index}) = 1;

HeartVoxels_cut = 1-inv;
HeartVoxels = HeartVoxels_cut;

figure(12)
subplot(1,3,1)
imagesc(squeeze(HeartVoxels(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxels(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxels(end/2,:,:)));
daspect([1 2 1]);


figure(13)
subplot(1,3,1)
imagesc(squeeze(HeartVoxels(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxels(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxels(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);


%% Remove Ventricles

fprintf('removing ventricles...\n');

for i = 1:size(HeartVoxels,3)
   HeartVoxels(:,:,i) = imfill(squeeze(HeartVoxels(:,:,i)),'holes');     
end

for i = 1:size(HeartVoxels,1)
   HeartVoxels(i,:,:) = imfill(squeeze(HeartVoxels(i,:,:)),'holes');     
end

for i = 1:size(HeartVoxels,2)
   HeartVoxels(:,i,:) = imfill(squeeze(HeartVoxels(:,i,:)),'holes');     
end


figure(14)
subplot(1,3,1)
imagesc(squeeze(HeartVoxels(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxels(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxels(end/2,:,:)));
daspect([1 2 1]);

figure(15)
subplot(1,3,1)
imagesc(squeeze(HeartVoxels(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxels(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxels(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);



%% Remove Remaining Small Loose Parts

for i = 1:size(HeartVoxels,3)
    cc = bwconncomp(squeeze(HeartVoxels(:,:,i)));
    % find largest connected component
    [~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
    %rebuild component
    temp = zeros(size(HeartVoxels,1),size(HeartVoxels,2));
    if (~isempty(index))
        temp(cc.PixelIdxList{index}) = 1;
    end
    HeartVoxels(:,:,i) = temp.*squeeze(HeartVoxels(:,:,i));
end

for i = 1:size(HeartVoxels,1)
    cc = bwconncomp(squeeze(HeartVoxels(i,:,:)));
    % find largest connected component
    [~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
    %rebuild component
    temp = zeros(size(HeartVoxels,2),size(HeartVoxels,3));
    if (~isempty(index))
        temp(cc.PixelIdxList{index}) = 1;
    end
    HeartVoxels(i,:,:) = temp.*squeeze(HeartVoxels(i,:,:));
end

for i = 1:size(HeartVoxels,2)
    cc = bwconncomp(squeeze(HeartVoxels(:,i,:)));
    % find largest connected component
    [~,index] = max(cellfun(@(x) numel(x),cc.PixelIdxList));
    %rebuild component
    temp = zeros(size(HeartVoxels,1),size(HeartVoxels,3));
    if (~isempty(index))
        temp(cc.PixelIdxList{index}) = 1;
    end
    HeartVoxels(:,i,:) = temp.*squeeze(HeartVoxels(:,i,:));
end

figure(16)
subplot(1,3,1)
imagesc(squeeze(HeartVoxels(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxels(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxels(end/2,:,:)));
daspect([1 2 1]);

figure(17)
subplot(1,3,1)
imagesc(squeeze(HeartVoxels(:,:,end/2)).*squeeze(images_TE1(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxels(:,end/2,:)).*squeeze(images_TE1(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxels(end/2,:,:)).*squeeze(images_TE1(end/2,:,:)));
daspect([1 2 1]);


%% Create Shell
fprintf('creating shell...\n');

shell_thick = 3;
HeartVoxelsShell = zeros(size(HeartVoxels));

se = strel('disk',round(shell_thick/0.2879),0);
for i = 1:size(HeartVoxels,3)
   HeartVoxelsShell(:,:,i) = imfill(round(imdilate(squeeze(HeartVoxels(:,:,i)),se)));     
end

se = strel('line',round(shell_thick/0.58),0);
for i = 1:size(HeartVoxels,2)
   HeartVoxelsShell(:,i,:) = imfill(round(imdilate(squeeze(HeartVoxelsShell(:,i,:)),se)));     
end

HeartVoxelsShell(HeartVoxels==1)=0;
HeartVoxelsShell(:,:,1:cut_slice)=0;

figure(18)
subplot(1,3,1)
imagesc(squeeze(HeartVoxelsShell(:,:,end/2)));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxelsShell(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxelsShell(end/2,:,:)));
daspect([1 2 1]);


%% Create stand
lxy = 25;
lz = 50;

maxs = squeeze(max(max(HeartVoxelsShell)));
increase = 1:size(HeartVoxelsShell,3);
z_arr = maxs.*increase';
z_apex = max(z_arr);

if z_apex>size(HeartVoxels,3)-3
    HeartVoxelsShell(:,:,end:end+5)=0;
end

stand = zeros(size(HeartVoxelsShell));
stand(round(end/2-20-lxy/0.2879):round(end/2-20+lxy/0.2879),round(end/2+20-lxy/0.2879):round(end/2+20+lxy/0.2879),round(z_apex-lz/0.58):z_apex+3) = 1;
stand(:,:,1:size(HeartVoxels,3)) = stand(:,:,1:size(HeartVoxels,3))-HeartVoxels;
stand(stand<0)=0;
HeartVoxelsShellStand = HeartVoxelsShell;
HeartVoxelsShellStand(stand==1)=1;

figure(19)
subplot(1,3,1)
imagesc(squeeze(HeartVoxelsShellStand(:,:,round(end/2))));
daspect([1 1 1]);
subplot(1,3,2)
imagesc(squeeze(HeartVoxelsShellStand(:,end/2,:)));
daspect([1 2 1]);
subplot(1,3,3)
imagesc(squeeze(HeartVoxelsShellStand(end/2,:,:)));
daspect([1 2 1]);

%% Imprint name

rgb = insertText(zeros(size(HeartVoxels,1), size(HeartVoxels,2)), 1/2*[size(HeartVoxels,1) size(HeartVoxels,2)], 'Test', ...
        'BoxColor',     'black', ...
        'BoxOpacity',   0, ...
        'TextColor',    'white', ...
        'FontSize',     50, ...
        'AnchorPoint',  'Center');

textmask_2d = round(rgb(:, :, 1));

textmask = zeros(size(HeartVoxelsShellStand));
textmask(:,:,z_apex+2) = flipud(textmask_2d);
textmask(:,:,z_apex+3) = flipud(textmask_2d);
textmask(:,:,z_apex+4) = flipud(textmask_2d);

HeartVoxelsShellStandText = HeartVoxelsShellStand;
HeartVoxelsShellStandText(textmask==1)=0;

%% Create Heart Surface

[x, y, z] = meshgrid(1:size(HeartVoxels,1),1:size(HeartVoxels,2),1:size(HeartVoxels,3));

Iso = isosurface(0.2879.*x,0.2879.*y,0.58.*z,HeartVoxels,0.9);
sIso = smoothpatch(Iso,1,2);

%% Set Targets

t1 = [315,191,58];
t2 = [308,190,81];
t3 = [286,184,98];

%%

targets = [t1;t2;t3];

line_mask = zeros(2.*[size(mask_final_cut,1) size(mask_final_cut,2) 2*size(mask_final_cut,3)]);
guide_mask = zeros(2.*[size(mask_final_cut,1) size(mask_final_cut,2) 2*size(mask_final_cut,3)]);
distances = zeros(size(targets,1),1);

for targ = 1:size(targets,1)
    
    x_targ = targets(targ,1);
    y_targ = targets(targ,2);
    z_targ = targets(targ,3);

    points = [y_targ*0.2879 x_targ*0.2879 z_targ*0.58];
    [distances(targ),surface_points] = point2trimesh(sIso, 'QueryPoints', points,'UseSubSurface',false);
    

    y_insert = round((surface_points(1))/0.2879);
    x_insert = round((surface_points(2))/0.2879);
    z_insert = round((surface_points(3))/0.58);


    y_insert2 = 2*y_insert;
    x_insert2 = 2*x_insert;
    z_insert2 = 2*z_insert;

    y_targ2 = 2*y_targ;
    x_targ2 = 2*x_targ;
    z_targ2 = 2*z_targ;

    steps = 1:400;
    steps = steps./max(steps);
    z_insert_temp = 2*z_insert2;
    z_targ_temp = 2*z_targ2;

    l1 = norm([x_insert2*0.2879/2, y_insert2*0.2879/2, z_insert_temp*0.29/2]-[x_targ2*0.2879/2, y_targ2*0.2879/2, z_targ_temp*0.29/2])+10;
    factor1 = l1/(norm([x_insert2*0.2879/2, y_insert2*0.2879/2, z_insert_temp*0.29/2]-[x_targ2*0.2879/2, y_targ2*0.2879/2, z_targ_temp*0.29/2]));

    distance = norm([x_insert2*0.2879/2, y_insert2*0.2879/2, z_insert_temp*0.29/2]-[x_targ2*0.2879/2, y_targ2*0.2879/2, z_targ_temp*0.29/2]);

    l2 = 2;
    factor2 = l2/(norm([x_insert2*0.2879/2, y_insert2*0.2879/2, z_insert_temp*0.29/2]-[x_targ2*0.2879/2, y_targ2*0.2879/2, z_targ_temp*0.29/2]));
    
    l3 = 3;
    factor3 = l3/(norm([x_insert2*0.2879/2, y_insert2*0.2879/2, z_insert_temp*0.29/2]-[x_targ2*0.2879/2, y_targ2*0.2879/2, z_targ_temp*0.29/2]));
    
    l4 = 20;
    factor4 = l4/(norm([x_insert2*0.2879/2, y_insert2*0.2879/2, z_insert_temp*0.29/2]-[x_targ2*0.2879/2, y_targ2*0.2879/2, z_targ_temp*0.29/2]));

    line_array1 = [x_targ2, y_targ2, z_targ_temp]+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp]).*factor1.*steps';
    line_array2 = [x_targ2, y_targ2, z_targ_temp]+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp])+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp])*factor2+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp]).*factor3.*steps';
    line_array3 = [x_targ2, y_targ2, z_targ_temp]+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp])+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp])*factor2+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp]).*factor3+([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp]).*factor4.*steps';

    traj_indices = ([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp])./norm([x_insert2, y_insert2, z_insert_temp]-[x_targ2, y_targ2, z_targ_temp]);

    orth = null(traj_indices(:).');
    orth1 = orth(:,1);
    orth2 = orth(:,2);

    for k = 1:0.5:360
        vec = rotVecAroundArbAxis(orth1',traj_indices,k);
        
        for i = 1:length(steps)    
            for j = -50:0.5:50               
                line_mask(round(line_array1(i,1)+vec(1)*j/50*4.9/2/0.2879*2),round(line_array1(i,2)+vec(2)*j/50*4.9/2/0.2879*2),round(line_array1(i,3)+vec(3)*j/50*4.9/2/0.29*2)) = 1;
                
                if (round(line_array3(i,1)+vec(1)*j/50*15/2/0.2879*2)<size(line_mask,1)&&round(line_array3(i,2)+vec(2)*j/50*15/2/0.2879*2)<size(line_mask,2)&&round(line_array3(i,3)+vec(3)*j/50*15/2/0.29*2)<size(line_mask,3))                                    
                    line_mask(round(line_array3(i,1)+vec(1)*j/50*15/2/0.2879*2),round(line_array3(i,2)+vec(2)*j/50*15/2/0.2879*2),round(line_array3(i,3)+vec(3)*j/50*15/2/0.29*2)) = 1;
                end
                
                if (round(line_array2(i,1)+vec(1)*j/50*15/2/0.2879*2)<size(guide_mask,1)&&round(line_array2(i,2)+vec(2)*j/50*15/2/0.2879*2)<size(guide_mask,2)&&round(line_array2(i,3)+vec(3)*j/50*15/2/0.29*2)<size(guide_mask,3))                                    
                    guide_mask(round(line_array2(i,1)+vec(1)*j/50*15/2/0.2879*2),round(line_array2(i,2)+vec(2)*j/50*15/2/0.2879*2),round(line_array2(i,3)+vec(3)*j/50*15/2/0.29*2)) = 1;
                end
            end
        end
    end
    
end

distances = abs(distances);

%% Boolean Combine/Subtract of Masks

line_mask_res = imresize3(line_mask,size(HeartVoxelsShellStandText));
guide_mask_res = imresize3(guide_mask,size(HeartVoxelsShellStandText));

final_mask = HeartVoxelsShellStandText;
final_mask(round(imresize3(guide_mask,size(HeartVoxelsShellStandText)))==1)=1;
final_mask(round(line_mask_res)==1)=0;

%% Create and Export STL

[x, y, z] = meshgrid(1:size(final_mask,1),1:size(final_mask,2),1:size(final_mask,3));

Iso = isosurface(0.2879.*x,0.2879.*y,0.58.*z,final_mask,0.9);

% center around origin
X = [min(Iso.vertices(:,1)) max(Iso.vertices(:,1))];
Y = [min(Iso.vertices(:,2)) max(Iso.vertices(:,2))];
Z = [min(Iso.vertices(:,3)) max(Iso.vertices(:,3))];

Iso.vertices(:,1) = Iso.vertices(:,1) - (X(1)+X(2))/2;
Iso.vertices(:,2) = Iso.vertices(:,2) - (Y(1)+Y(2))/2;
Iso.vertices(:,3) = Iso.vertices(:,3) - (Z(1)+Z(2))/2;


fprintf('Smoothing 3D Model...\n');
% smoothing, 3 Iterations
sIso = smoothpatch(Iso,1,3);

figure(20)
p = patch(sIso);
isonormals(0.2879.*x,0.2879.*y,0.58.*z,final_mask, p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3)
camlight; lighting phong

TR = triangulation(sIso.faces,sIso.vertices);

fprintf('Saving 3D Model...\n');
stlwrite(TR,'Heart.stl');


