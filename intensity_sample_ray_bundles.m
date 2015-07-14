function [] = intensity_sample_ray_bundles(im_filename, label, r_center, c_center)

%% setup

batch_mode               = 1;        % in batch mode, figures are invisible and no output is displayed
im.raw                   = imread(im_filename);  % raw image filename
r_dim                    = size(im.raw, 1);  % row dimension of the raw image
c_dim                    = size(im.raw, 2);  % column dimension of the raw image
smooth_num               = 100;      % number of times to smooth the gray image
smooth_mask              = [5 5];    % size of mask to smooth the gray image
smooth_std               = 5.0;      % standard deviation of the Gaussian function used to smooth the gray image
delta_theta_bundle       = 2*pi;     % angle of bundle: 2*pi (full circle => one bundle), pi/4 (1/8 circle => 8 bundles), etc.
delta_theta              = pi/40;    % angle between each extended sample ray: pi/40 => 80 rays per full circle, etc.
delta_rho                = 1;        % radial distance between intensity level samples (number of pixels)
delta_n                  = 0;        % radius of neighborhood over which to average intensity levels at each sample
radius_mean_min_window   = 1000;     % maximum radial distance over which to compute the global minimum mean for determining the radius of each bundle
radius_median_min_window = 1000;     % maximum radial distance over which to compute the global minimum median for determining the radius of each bundle
max_rho_sample           = ceil(sqrt(r_dim^2+c_dim^2)/delta_rho);  % maximum number of samples that each sample ray can take (for preallocating arrays)
max_theta_sample         = ceil(2*pi/delta_theta);                 % maximum number of sample rays that can fit in a circle (for preallocating arrays, setting index limits)
max_theta_bundle         = ceil(2*pi/delta_theta_bundle);          % maximum number of bundles that can fit in a circle (for setting index limits)
bundle_samples           = ceil(delta_theta_bundle/delta_theta);   % number of sample rays that can fit in a bundle (for determining ranges)
center_num               = numel(r_center);  % number of centers to process

% create new results directory
path = strcat('/tmp/analysis/gradient/', label);
mkdir(path);

% open a data file
data_file = strcat(path, '/', 'data.txt');
fid = fopen(data_file, 'w');

fprintf(fid, 'number of bundles per circle = %f\n', max_theta_bundle);
fprintf(fid, 'number of sample rays per circle = %f\n', max_theta_sample);
fprintf(fid, '\n');

% convert RGB image to gray image
im.gray = rgb2gray(im.raw);

% smooth gray image
G = fspecial('gaussian', smooth_mask, smooth_std);
im.smooth = im.gray;
for i = 1 : smooth_num
    im.smooth = imfilter(im.smooth, G, 'replicate', 'conv');
end

% for each center, measure intensity levels using bundles of sample rays, and plot the quantitative results for each bundle
for center_iter = 1 : center_num
    
    fprintf(fid, 'center %d\n', center_iter);
    
    % setup
    r_init       = r_center(center_iter);
    c_init       = c_center(center_iter);
    measurements = 256 * ones(max_rho_sample, max_theta_sample);
    
    % measure
    for theta_iter = 1 : max_theta_sample
        r        = r_init;
        c        = c_init;
        theta    = (theta_iter - 1)*delta_theta;
        rho      = 0;
        rho_iter = 1;
        while r >= 1 && r <= r_dim && c >= 1 && c <= c_dim;
            neighbor_val = [];
            neighbor_num = 1;
            for rr = -delta_n : delta_n
                for cc = -delta_n : delta_n
                    if r+rr >= 1 && r+rr <= r_dim && c+cc >=1 && c+cc <= c_dim
                        neighbor_val(neighbor_num) = im.smooth(r+rr,c+cc);
                        neighbor_num = neighbor_num + 1;
                    end
                end
            end
            val = mean(neighbor_val);
            measurements(rho_iter,theta_iter) = val;
            if ~batch_mode
                fprintf('%d/%d: (%d,%f) = (%d,%d): %d\n', center_iter, center_num, rho, theta/pi, r, c, val);
            end
            rho = rho + delta_rho;
            [x, y] = pol2cart(theta, rho);
            r = r_init - floor(y);
            c = c_init + floor(x);
            rho_iter = rho_iter + 1;
        end
    end
    
    % compute statistics and plot
    for bundle_iter = 1 : max_theta_bundle
        
        fprintf(fid, '\tbundle %d\n', bundle_iter);
        
        measurements_bundle_range      = (bundle_iter - 1)*bundle_samples + 1 : bundle_iter*bundle_samples;
        measurements_mean              = safe_stats(@mean, measurements(:,measurements_bundle_range)); % mean(measurements(:,measurements_bundle_range), 2);
        measurements_std               = safe_stats(@std, measurements(:,measurements_bundle_range)); % std(measurements(:,measurements_bundle_range), 0, 2);
        measurements_cv                = measurements_std ./ measurements_mean;
        measurements_median            = safe_stats(@median, measurements(:,measurements_bundle_range)); % median(measurements(:,measurements_bundle_range), 2);
        measurements_median_min_locs   = find_mins(measurements_median, radius_median_min_window);
        measurements_radius            = measurements_median_min_locs(1) * delta_rho;
        measurements_range             = 1:measurements_radius;
        radii(center_iter,bundle_iter) = measurements_radius;

        fprintf(fid, '\t\tradius = %d\n', measurements_radius);
        
        if batch_mode
            h = figure('Visible', 'off');
        else
            h = figure(bundle_iter);
        end
        clf;
        
        subplot_width = 3;

        % plot all trajectories, overlaid with mean and median
        subplot(1,subplot_width,1);
        hold on;
        safe_plot(measurements(measurements_range,:));
        plot(measurements_mean(measurements_range), 'b');
        plot(measurements_median(measurements_range), 'r');
        hold off;
        title(['radius=', num2str(measurements_radius)]);
        
        % plot the mean +/- std
        fprintf(fid, '\t\tmean statistics\n');
        subplot(1,subplot_width,2);
        plot_time_series(measurements_mean', measurements_std', [0 0 1], 1, measurements_radius);
        hold on;
        [seg_i, seg_j] = segmented_least_squares(measurements_range, measurements_mean(measurements_range)', 200);
        num_segs = numel(seg_i);
        seg_a = zeros(1,num_segs);
        seg_b = zeros(1,num_segs);
        seg_e = zeros(1,num_segs);
        for s = 1 : num_segs
            fprintf(fid, '\t\t\tsegment %d\n', s);
            [seg_a(s), seg_b(s)] = least_squares_fit(measurements_range, measurements_mean(measurements_range)', seg_i(s), seg_j(s));
            seg_e(s) = least_squares_error(measurements_range, measurements_mean(measurements_range)', seg_i(s), seg_j(s));
            plot(measurements_range(seg_i(s):seg_j(s)), seg_a(s) .* measurements_range(seg_i(s):seg_j(s)) + seg_b(s), 'k');
            fprintf(fid, '\t\t\t\tbeg = %d\n', seg_i(s));
            fprintf(fid, '\t\t\t\tend = %d\n', seg_j(s));
            fprintf(fid, '\t\t\t\tlen = %d\n', seg_j(s) - seg_i(s) + 1);
            fprintf(fid, '\t\t\t\tslo = %f\n', seg_a(s));
            fprintf(fid, '\t\t\t\terr = %f\n', seg_e(s));
        end
        hold off;
        title_foo = {};
        for s = 1 : num_segs
            title_foo{numel(title_foo)+1} = ['l=', num2str(seg_j(s) - seg_i(s) + 1), ', s=', sprintf('%0.2f', seg_a(s)), ', e=', sprintf('%0.2f', seg_e(s))];
        end
        title(title_foo);
        
        % plot the median +/- std
        fprintf(fid, '\t\tmedian statistics\n');
        subplot(1,subplot_width,3);
        plot_time_series(measurements_median', measurements_std', [1 0 0], 1, measurements_radius);
        hold on;
        [seg_i, seg_j] = segmented_least_squares(measurements_range, measurements_median(measurements_range)', 200);
        num_segs = numel(seg_i);
        seg_a = zeros(1,num_segs);
        seg_b = zeros(1,num_segs);
        seg_e = zeros(1,num_segs);
        for s = 1 : numel(seg_i)
            fprintf(fid, '\t\t\tsegment %d\n', s);
            [seg_a(s), seg_b(s)] = least_squares_fit(measurements_range, measurements_median(measurements_range)', seg_i(s), seg_j(s));
            seg_e(s) = least_squares_error(measurements_range, measurements_mean(measurements_range)', seg_i(s), seg_j(s));
            plot(measurements_range(seg_i(s):seg_j(s)), seg_a(s) .* measurements_range(seg_i(s):seg_j(s)) + seg_b(s), 'k');
            fprintf(fid, '\t\t\t\tbeg = %d\n', seg_i(s));
            fprintf(fid, '\t\t\t\tend = %d\n', seg_j(s));
            fprintf(fid, '\t\t\t\tlen = %d\n', seg_j(s) - seg_i(s) + 1);
            fprintf(fid, '\t\t\t\tslo = %f\n', seg_a(s));
            fprintf(fid, '\t\t\t\terr = %f\n', seg_e(s));
        end
        hold off;
        title_foo = {};
        for s = 1 : num_segs
            title_foo{numel(title_foo)+1} = ['l=', num2str(seg_j(s) - seg_i(s) + 1), ', s=', sprintf('%0.2f', seg_a(s)), ', e=', sprintf('%0.2f', seg_e(s))];
        end
        title(title_foo);
        
        fig_file = strcat(path, '/', 'blur_radius_', num2str(center_iter), '_bundle_', num2str(bundle_iter), '.pdf');
        saveas(gcf, fig_file, 'pdf');
        
    end
    
end

% plot the smoothed image, overlaid with centers, bundles, and bundle numbers
if batch_mode
    h = figure('Visible', 'off');
else
    h = figure(bundle_iter);
end
clf;
imshow(im.smooth);
hold on;
for center_iter = 1 : numel(r_center)
    beg_x = [];
    beg_y = [];
    end_x = [];
    end_y = [];
    plot(c_center(center_iter), r_center(center_iter), 'ro');
    text(c_center(center_iter) + 30, r_center(center_iter), num2str(center_iter), 'Color', 'r');
    for bundle_iter = 1 : max_theta_bundle
        theta_range = (bundle_iter - 1)*delta_theta_bundle : 0.01 : bundle_iter*delta_theta_bundle;
        [perim_x, perim_y] = pol2cart(theta_range, radii(center_iter,bundle_iter));
        plot(c_center(center_iter) + perim_x, r_center(center_iter) - perim_y, 'r');
        [label_x, label_y] = pol2cart(theta_range(floor(numel(theta_range)/2)), radii(center_iter,bundle_iter) + 30);
        text(c_center(center_iter) + label_x, r_center(center_iter) - label_y, num2str(bundle_iter), 'Color', 'r');
        beg_x(bundle_iter) = c_center(center_iter) + perim_x(1);
        beg_y(bundle_iter) = r_center(center_iter) - perim_y(1);
        end_x(bundle_iter) = c_center(center_iter) + perim_x(numel(perim_x));
        end_y(bundle_iter) = r_center(center_iter) - perim_y(numel(perim_y));
    end
    for bundle_iter = 1 : max_theta_bundle - 1
        plot([end_x(bundle_iter) beg_x(bundle_iter+1)], [end_y(bundle_iter) beg_y(bundle_iter+1)], 'r');
    end
    plot([end_x(bundle_iter+1) beg_x(1)], [end_y(bundle_iter+1) beg_y(1)], 'r');
end
hold off;
fig_file = strcat(path, '/', 'blur_radii.pdf');
saveas(gcf, fig_file, 'pdf');

% close the data file
fclose(fid);

return;

    %% plot a time series average +/- standard deviation (average curve surrounded by +/- gray patches)
    %  a: average            (m time series x n ticks)
    %  s: standard deviation (m time series x n ticks)
    %  c: color map          (m time series x 3 {r,g,b})
    %  b: begin time         (integer)
    %  e: end time           (integer)
    function [] = plot_time_series(a, s, c, b, e)
        domain = b : e;
        gray = [0.9 0.9 0.9];
        hold on;
        for t = 1 : size(a,1)
            patch([domain fliplr(domain)], [a(t,domain) - s(t,domain), fliplr(a(t,domain) + s(t,domain))], gray, 'LineStyle', 'none');
            plot(domain, a(t,domain), 'color', c(t,:));
        end
        hold off;
    end  % function plot_time_series

    function [l] = find_mins(a, w)
        [m, l] = min(a(1:w));
    end

    function [F] = safe_stats(f, M)
        num_r = size(M,1);
        num_c = size(M,2);
        F = [];
        for r = 1 : num_r
            safe_set = [];
            for c = 1 : num_c
                if M(r,c) ~= 256
                    safe_set(numel(safe_set)+1) = M(r,c);
                end
            end
            F(r,1) = f(safe_set);
        end
    end

    function [] = safe_plot(M)
        num_r = size(M,1);
        num_c = size(M,2);
        for c = 1 : num_c
            for r = 1 : num_r
                if M(r,c) == 256
                    break;
                end
            end
            safe_r_idx = r - 1;
            plot(1:safe_r_idx, M(1:safe_r_idx, c), 'Color', [0.9 0.9 0.9]);
        end
    end

end
