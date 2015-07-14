function [] = quad_tree(im_filename)

%% load image
im.raw  = imread(im_filename);
im.gray = rgb2gray(im.raw);
r_dim = size(im.raw, 1);
c_dim = size(im.raw, 2);

%% smooth image
smooth_num = 100;
G = fspecial('gaussian', [5 5], 5.0);
im.smooth = im.gray;
for i = 1 : smooth_num
    im.smooth = imfilter(im.smooth, G, 'replicate', 'same', 'conv');
end

%% quad tree
r_range = 1:r_dim;
c_range = 1:c_dim;
I = im.smooth(r_range,c_range);

fh=figure;
imshow(I);
hold on;
dissect(I, 1, 1, 0, @(x) std(x)/mean(x), 0.02, 1);
hold off;

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
        end     
    end

end
