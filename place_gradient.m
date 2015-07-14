function [] = place_gradient(he_filename, pimo_filename, label, pimo_center_rows, pimo_center_cols, he_center_rows, he_center_cols)

%% parameters
theta_inc   = 0.01;
theta_range = 0 : theta_inc : 2*pi + theta_inc;

%% load images
base_path = '/tmp/images';
full_path = strcat(base_path, '/', label);
he_file   = strcat(full_path, '/', he_filename);
pimo_file = strcat(full_path, '/', pimo_filename);
im.he_raw    = imread(he_file);  % H&E image
im.pimo_raw  = imread(pimo_file);  % anti-pimo image
im.pimo_gray = rgb2gray(im.pimo_raw);
r_dim = size(im.pimo_raw, 1);
c_dim = size(im.pimo_raw, 2);

% row beg and end index values into the gradient matrix below
gradient_idx = [1 3  ;
                4 8  ;
                9 11 ];

% gradient matrix: delta r, slope, initial intensity
gradient     = [ 70 -0.38 204 ;
                111 -0.00 174 ;
                123 -0.15 173 ;
                 51 -0.44 200 ;
                 88 -0.05 178 ;
                110 -0.09 177 ;
                 77 -0.24 169 ;
                 36  0.05 151 ;
                 68 -0.28 197 ;
                108 -0.06 178 ;
                144 -0.15 175 ];

%% smooth image
smooth_num = 100;
G = fspecial('gaussian', [5 5], 5.0);
im.pimo_smooth = im.pimo_gray;
for i = 1 : smooth_num
    im.pimo_smooth = imfilter(im.pimo_smooth, G, 'replicate', 'conv');
end
im.pimo_smooth_flip = im.pimo_smooth(1920:-1:1,:);

%% set up figure
h1 = figure(1);
clf;
hold on;

%% plot centers and perimeters
for center_idx = 1 : numel(radii)
    plot( pimo_center_cols(center_idx), pimo_center_rows(center_idx), 'ro' );
    [perim_x, perim_y] = pol2cart(theta_range, radii(center_idx));
    plot(pimo_center_cols(center_idx) + perim_x, pimo_center_rows(center_idx) - perim_y, 'r');
end

%% plot concentric rings based on gradient intensity values
for center_idx = 1 : numel(radii)
    beg_grad_idx = gradient_idx(center_idx, 1);
    end_grad_idx = gradient_idx(center_idx, 2);
    r_cursor = 0;
    for grad_idx = beg_grad_idx : end_grad_idx
        init_r   = r_cursor;
        delta_r  = gradient(grad_idx, 1);
        slope    = gradient(grad_idx, 2);
        init_i   = gradient(grad_idx, 3);
        i_cursor = init_i;
        for r_idx = r_cursor : init_r + delta_r
            r_cursor = r_cursor + 1;
            [perim_x, perim_y] = pol2cart(theta_range, r_cursor);
            i_cursor_norm = i_cursor / 255;
            plot(pimo_center_cols(center_idx) + perim_x, pimo_center_rows(center_idx) - perim_y, 'Color', [i_cursor_norm, i_cursor_norm, i_cursor_norm]);
            i_cursor = i_cursor + slope;
        end
    end
end

%% show raw pimo image (at low opacity)
imshow(im.pimo_raw);
alpha(0.5);

% enlarge and save the figure
set(gcf, 'PaperPosition', [0 0 20 20])    % can be bigger than screen 
set(gcf, 'PaperSize', [20 20])    % Same, but for PDF output
pimo_gradient_overlay_file = strcat(full_path, '/', 'pimo_gradient_overlay.pdf');
print(gcf, pimo_gradient_overlay_file, '-dpdf', '-r300' );

%% set up figure
h2 = figure(2);
clf;
hold on;

%% plot centers and perimeters
for center_idx = 1 : numel(radii)
    plot( he_center_cols(center_idx), he_center_rows(center_idx), 'ro' );
    [perim_x, perim_y] = pol2cart(theta_range, radii(center_idx));
    plot(he_center_cols(center_idx) + perim_x, he_center_rows(center_idx) - perim_y, 'r');
end

%% plot concentric rings based on gradient intensity values
for center_idx = 1 : numel(radii)
    beg_grad_idx = gradient_idx(center_idx, 1);
    end_grad_idx = gradient_idx(center_idx, 2);
    r_cursor = 0;
    for grad_idx = beg_grad_idx : end_grad_idx
        init_r   = r_cursor;
        delta_r  = gradient(grad_idx, 1);
        slope    = gradient(grad_idx, 2);
        init_i   = gradient(grad_idx, 3);
        i_cursor = init_i;
        for r_idx = r_cursor : init_r + delta_r
            r_cursor = r_cursor + 1;
            [perim_x, perim_y] = pol2cart(theta_range, r_cursor);
            i_cursor_norm = i_cursor / 255;
            plot(he_center_cols(center_idx) + perim_x, he_center_rows(center_idx) - perim_y, 'Color', [i_cursor_norm, i_cursor_norm, i_cursor_norm]);
            i_cursor = i_cursor + slope;
        end
    end
end

%% show raw H&E image (at low opacity)
imshow(im.he_raw);
alpha(0.5);

% enlarge and save the figure
set(gcf, 'PaperPosition', [0 0 20 20])    % can be bigger than screen 
set(gcf, 'PaperSize', [20 20])    % Same, but for PDF output
he_gradient_overlay_file = strcat(full_path, '/', 'he_gradient_overlay.pdf');
print(gcf, he_gradient_overlay_file, '-dpdf', '-r300' );

%% plot H&E -> pimo center vectors
h3 = figure(3);
clf;
imshow(im.pimo_raw);
hold on;
for center_idx = 1 : numel(radii)
    plot(he_center_cols(center_idx), he_center_rows(center_idx), 'ko');
    plot(pimo_center_cols(center_idx), pimo_center_rows(center_idx), 'ko');
    quiver(he_center_cols(center_idx), he_center_rows(center_idx), pimo_center_cols(center_idx) - he_center_cols(center_idx), pimo_center_rows(center_idx) - he_center_rows(center_idx), 0.98, 'b');

    he_y_buffer   = 25;
    pimo_y_buffer = 25;
    if center_idx == 3
        he_y_buffer   = -he_y_buffer;
        pimo_y_buffer = -pimo_y_buffer;
    end
    x_buffer = 70;
    
    text(he_center_cols(center_idx), he_center_rows(center_idx) + he_y_buffer, strcat('H_{', num2str(center_idx), '}'), 'Color', 'b');
    text(pimo_center_cols(center_idx), pimo_center_rows(center_idx) + pimo_y_buffer, strcat('P_{', num2str(center_idx), '}'), 'Color', 'b');

    dy(center_idx)    = he_center_rows(center_idx) - pimo_center_rows(center_idx);
    dx(center_idx)    = he_center_cols(center_idx) - pimo_center_cols(center_idx);
    r(center_idx)     = sqrt(dy(center_idx)^2 + dx(center_idx)^2);
    theta(center_idx) = asin(dy(center_idx) / r(center_idx));
    fprintf('vec %d: (%d,%d) -> (%d,%d), dx = %d, dy = %d, r = %f, theta = %f\n', center_idx, he_center_cols(center_idx), he_center_rows(center_idx), pimo_center_cols(center_idx), pimo_center_rows(center_idx), dx(center_idx), dy(center_idx), r(center_idx), theta(center_idx));

    range = he_center_cols(center_idx) : 10 : he_center_cols(center_idx) + r(center_idx);
    plot(range, he_center_rows(center_idx), 'b.', 'MarkerSize', 1);
    
    text(he_center_cols(center_idx) + x_buffer, he_center_rows(center_idx) - (he_y_buffer * 2/3), sprintf('%-1.3f rad', theta(center_idx)), 'Color', 'b');
    
    rot                = @(alpha,v) [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)] * v;
    text_init_pos      = [x_buffer ; -(he_y_buffer * 2/3)];
    text_rot_init_pos  = rot(theta(center_idx), text_init_pos);
    ht = text(he_center_cols(center_idx) + text_rot_init_pos(1), he_center_rows(center_idx) - text_rot_init_pos(2), sprintf('%-1.3f pix', r(center_idx)), 'Color', 'b');
    set(ht, 'rotation', theta(center_idx) * (360/(2*pi)));
end

% derive rotation angle theta as a function of x (rows)
P_theta = polyfit(he_center_cols, theta, 2);
F_theta = @(x) P_theta(1).*x.^2 + P_theta(2).*x + P_theta(3);

% derive a redius length as a function of y (cols)
P_radius = polyfit(he_center_rows, r, 2);
F_radius = @(y) P_radius(1).*y.^2 + P_radius(2).*y + P_radius(3);

% place vector field onto the image
for y_idx = 0 : 150 : 1900
    radius = F_radius(y_idx);
    for x_idx = 0 : 300 : 2500
        theta_idx = F_theta(x_idx);
        dy = -radius * sin(theta_idx);
        dx =  radius * cos(theta_idx);
        quiver(x_idx, y_idx, dx, dy, 1.0, 'r');
    end
end

% enlarge and save the figure
set(gcf, 'PaperPosition', [0 0 20 20])    % can be bigger than screen 
set(gcf, 'PaperSize', [20 20])    % Same, but for PDF output
he_to_pimo_vector_field_file = strcat(full_path, '/', 'he_to_pimo_vector_field.pdf');
print(gcf, he_to_pimo_vector_field_file, '-dpdf', '-r300' );

end  % end function place_gradient
