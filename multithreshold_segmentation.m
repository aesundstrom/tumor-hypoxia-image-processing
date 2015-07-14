function [] = multithreshold_segmentation(im_filename, label)

%% setup

batch_mode  = 1;  % in batch mode, figures are invisible and no output is displayed
im.raw      = imread(im_filename);  % raw image
r_dim       = size(im.raw, 1);  % row dimension of the raw image
c_dim       = size(im.raw, 2);  % column dimension of the raw image
smooth_num  = 100;    % number of times to smooth the gray image
smooth_mask = [5 5];  % size of mask to smooth the gray image
smooth_std  = 5.0;    % standard deviation of the Gaussian function used to smooth the gray image

%% preprocessing

% create new results directory
base_path = '/tmp/analysis/segment';
full_path = strcat(base_path, '/', label);
if ~exist(base_path, 'dir')
    mkdir(base_path);
end
if ~exist(full_path, 'dir')
    mkdir(full_path);
end

% open global and local data files
glo_data_file = strcat(base_path, '/', 'data.txt');
glo_fid       = fopen(glo_data_file, 'a');
loc_data_file = strcat(full_path, '/', 'data.txt');
loc_fid       = fopen(loc_data_file, 'w');

% convert RGB image to gray image
im.gry = rgb2gray(im.raw);

% smooth gray image
G = fspecial('gaussian', smooth_mask, smooth_std);
im.smo = im.gry;
for i = 1 : smooth_num
    im.smo = imfilter(im.smo, G, 'replicate', 'conv');
end

%% segment gray image

thresh = multithresh(im.gry, 2);
im.seg = imquantize(im.gry, thresh);
im.hyp = im.seg == 1;
im.via = im.seg == 2;
im.nec = im.seg == 3;
num_pix = r_dim * c_dim;
num_hyp = sum(sum(im.hyp));
num_via = sum(sum(im.via));
num_nec = sum(sum(im.nec));
fprintf(glo_fid, '%s\t%d\t%d\t%d\t%d\t%d\t%d\t', label, num_pix, thresh(1), thresh(2), num_hyp, num_via, num_nec);
fprintf(loc_fid, '%s\t%d\t%d\t%d\t%d\t%d\t%d\t', label, num_pix, thresh(1), thresh(2), num_hyp, num_via, num_nec);

%% plot gray and segmented-gray images side-by-side

if batch_mode
    h = figure('Visible', 'off');
else
    h = figure;
end
im.rgb = label2rgb(im.seg);
imshowpair(im.gry, im.rgb, 'montage');
fig_file = strcat(full_path, '/', 'segmented_gray.pdf');
saveas(gcf, fig_file, 'pdf');

%% segment smooth image

thresh = multithresh(im.smo, 2);
im.seg = imquantize(im.smo, thresh);
im.hyp = im.seg == 1;
im.via = im.seg == 2;
im.nec = im.seg == 3;
num_pix = r_dim * c_dim;
num_hyp = sum(sum(im.hyp));
num_via = sum(sum(im.via));
num_nec = sum(sum(im.nec));
fprintf(glo_fid, '%d\t%d\t%d\t%d\t%d\n', thresh(1), thresh(2), num_hyp, num_via, num_nec);
fprintf(loc_fid, '%d\t%d\t%d\t%d\t%d\n', thresh(1), thresh(2), num_hyp, num_via, num_nec);

%% plot smooth and segmented-smooth images side-by-side

if batch_mode
    h = figure('Visible', 'off');
else
    h = figure;
end
im.rgb = label2rgb(im.seg);
imshowpair(im.smo, im.rgb, 'montage');
fig_file = strcat(full_path, '/', 'segmented_smooth.pdf');
saveas(gcf, fig_file, 'pdf');

% close the data files
fclose(glo_fid);
fclose(loc_fid);

return;

end
