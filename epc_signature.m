function [] = epc_signature(im_filename, label)

%% setup

batch_mode = 1;      % in batch mode, figures are invisible and no output is displayed
seg_cost   = 50000;  % cost to add a segment in the segmented least squares algorithm 
im.raw     = imread(im_filename);  % raw image
r_dim      = size(im.raw, 1);  % row dimension of the raw image
c_dim      = size(im.raw, 2);  % column dimension of the raw image
gray_lim   = 255;         % maximum intensity level in an 8-bit image
x          = 1:gray_lim;  % the intensity range in an 8-bit image

%% preprocessing

% create new results directory
base_path = '/tmp/analysis/epc';
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

%% compute EPC curve for gray image

fprintf(glo_fid, '%s\t', label);
fprintf(loc_fid, '%s\t', label);
chi = zeros(1,gray_lim);
for i = 1 : gray_lim
    im.bin = im2bw(im.gry, i/gray_lim);
    [chi(i), l] = imEuler2d(im.bin);
    fprintf(glo_fid, '%d\t', chi(i));
    fprintf(loc_fid, '%d\t', chi(i));
end

%% compute segmented least-squares fit
[seg_i, seg_j] = segmented_least_squares(x, chi, seg_cost);
num_segs = numel(seg_i);
seg_a = zeros(1,num_segs);
seg_b = zeros(1,num_segs);
seg_e = zeros(1,num_segs);
compression_factor = gray_lim/(2*num_segs);
fprintf(glo_fid, '%d\t%f\t', num_segs, compression_factor);
fprintf(loc_fid, '%d\t%f\t', num_segs, compression_factor);
for s = 1 : num_segs
    [seg_a(s), seg_b(s)] = least_squares_fit(x, chi, seg_i(s), seg_j(s));
    seg_e(s) = least_squares_error(x, chi, seg_i(s), seg_j(s));
end
tot_err  = sum(seg_e);
norm_err = tot_err / gray_lim;
fprintf(glo_fid, '%f\t%f\t', tot_err, norm_err);
fprintf(loc_fid, '%f\t%f\t', tot_err, norm_err);
for s = 1 : num_segs
    fprintf(glo_fid, '%f\t%f', seg_a(s), seg_b(s));
    fprintf(loc_fid, '%f\t%f', seg_a(s), seg_b(s));
    if s ~= num_segs
        fprintf(glo_fid, '\t');
        fprintf(loc_fid, '\t');
    end
end
fprintf(glo_fid, '\n');
fprintf(loc_fid, '\n');

%% plot EPC curve for gray image

if batch_mode
    h = figure('Visible', 'off');
else
    h = figure;
end
clf;
hold on;
plot(x, chi, 'b');
for s = 1 : numel(seg_i)
    plot(x(seg_i(s):seg_j(s)), seg_a(s) .* x(seg_i(s):seg_j(s)) + seg_b(s), 'r');
end
hold off;
fig_file = strcat(full_path, '/', 'epc_signature.pdf');
saveas(gcf, fig_file, 'pdf');

% close the data files
fclose(glo_fid);
fclose(loc_fid);

return;

end
