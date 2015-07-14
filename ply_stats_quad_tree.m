function [] = ply_stats_quad_tree(im_filename, label)

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
base_path = '/tmp/analysis/quadtree';
full_path = strcat(base_path, '/', label);
if ~exist(base_path, 'dir')
    mkdir(base_path);
end
if ~exist(full_path, 'dir')
    mkdir(full_path);
end

% open global and local data files
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

%% quad tree
r_range = 1:r_dim;
c_range = 1:c_dim;
I = im.smo(r_range,c_range);

%% ply-by-ply quadtree
%    ply 1 is whole image
%    ply 2 is 1/4 of whole image
%    ply q is (1/4)^(q-1) of whole image
%    for each trial value of CV in [0.01 : 0.01 : 0.10] report statistics for each ply's set of windows
ply.n      = zeros(1,10);
ply.sum    = zeros(1,10);
ply.mean   = zeros(1,10);
ply.median = zeros(1,10);
ply.std    = zeros(1,10);
ply.cv     = zeros(1,10);
ply.hist   = zeros(12,10);
for t_iter = 1 : 10
    t = t_iter * 0.01;
    fprintf('Processing t=%f\n', t);
    ply.data = [];
    dissect(I, 1, 1, 0, @(x) std(x)/mean(x), t, 0);
    ply.n(t_iter) = numel(ply.data);
    ply.sum(t_iter) = sum(ply.data);
    ply.mean(t_iter) = mean(ply.data);
    ply.median(t_iter) = median(ply.data);
    ply.std(t_iter) = std(ply.data);
    ply.cv(t_iter) = ply.std / ply.mean;
    for p_iter = 1 : 12
        ply.hist(t_iter,p_iter) = sum(ply.data == p_iter);
    end
end

%% write data

head = strcat(label, '\t');

for t_iter = 1 : 10
    glo_data_file = strcat(base_path, '/', 'data_', num2str(t_iter), '.txt');
    glo_fid = fopen(glo_data_file, 'a');
    fprintf(glo_fid, head);
    fprintf(glo_fid, '%d\t', ply.hist(t_iter,:));
    fprintf(glo_fid, '\n');
    fclose(glo_fid);
end

fprintf(loc_fid, head);
fprintf(loc_fid, '%d\t', ply.n);
fprintf(loc_fid, '\n');
fprintf(loc_fid, head);
fprintf(loc_fid, '%d\t', ply.sum);
fprintf(loc_fid, '\n');
fprintf(loc_fid, head);
fprintf(loc_fid, '%f\t', ply.mean);
fprintf(loc_fid, '\n');
fprintf(loc_fid, head);
fprintf(loc_fid, '%d\t', ply.median);
fprintf(loc_fid, '\n');
fprintf(loc_fid, head);
fprintf(loc_fid, '%f\t', ply.std);
fprintf(loc_fid, '\n');
fprintf(loc_fid, head);
fprintf(loc_fid, '%f\t', ply.cv);
fprintf(loc_fid, '\n');
for t_iter = 1 : 10
    fprintf(loc_fid, head);
    fprintf(loc_fid, '%d\t', ply.hist(t_iter,:));
    fprintf(loc_fid, '\n');
end

% close the local data file
fclose(loc_fid);

return;

    function [] = draw_cross(r, c, rs, cs)
        hrs = floor(rs/2);
        hcs = floor(cs/2);
        plot(c+hcs,  r:r+rs, 'r');  % NS
        plot(c:c+cs, r+hrs,  'r');  % EW
    end

    %% recursively dissect a rectangle
    % I: 2-D image matrix
    % r: row coordinate of NW corner of rectangle
    % c: col coordinate of NW corner of rectangle
    % d: depth of recursion (0-based ply)
    % f: function handle for rectangle pixel evaluation function
    % t: threshold value in excess which will trigger further dissection
    % x: draw cross predicate
    function [] = dissect(I, r, c, d, f, t, x)
        if 2^d > min(size(I))
            ply.data(numel(ply.data)+1) = d;
            return;
        end
        rs = floor(size(I,1) / 2^d);
        cs = floor(size(I,2) / 2^d);
        r_range = r : r+rs-1;
        c_range = c : c+cs-1;
        W = I(r_range,c_range);
        v = double(reshape(W, 1, rs*cs));
        if f(v) > t
            if x
                draw_cross(r, c, rs, cs);
            end
            hrs = floor(rs/2);
            hcs = floor(cs/2);
            dissect(I, r,     c,     d+1, f, t, x);
            dissect(I, r+hrs, c,     d+1, f, t, x);
            dissect(I, r,     c+hcs, d+1, f, t, x);
            dissect(I, r+hrs, c+hcs, d+1, f, t, x);
        else
            ply.data(numel(ply.data)+1) = d;
        end        
    end

end
