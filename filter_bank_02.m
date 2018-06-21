clear all;
x = imread('cameraman.tif');

% analysis part
%
% horizontal filtering with H~(z^-1) and G~(z^-1)
x_h_H_tilde = zeros(size(x));
x_h_G_tilde = zeros(size(x));
for i=1:size(x,1)
    x_h_H_tilde(i,:) = H_tilde(x(i,:));
    x_h_G_tilde(i,:) = G_tilde(x(i,:));
end

% show image
img01 = uint8(x_h_H_tilde);
img02 = uint8(x_h_G_tilde);

% horizontal downsampling with a factor of 2
x_h_H_down = zeros(size(x,1),ceil(size(x,2)/2));
x_h_G_down = zeros(size(x,1),ceil(size(x,2)/2));
for i=1:size(x,1)
    x_h_H_down(i,:) = downsample(x_h_H_tilde(i,:),2);
    x_h_G_down(i,:) = downsample(x_h_G_tilde(i,:),2);
end

% vertical filtering
x_h_H_down_H_tilde = zeros(size(x_h_H_down));
x_h_H_down_G_tilde = zeros(size(x_h_H_down));
for i=1:size( x_h_H_down,2)
    x_h_H_down_H_tilde(:,i) = H_tilde(x_h_H_down(:,i)')';
    x_h_H_down_G_tilde(:,i) = G_tilde(x_h_H_down(:,i)')';
end

x_h_G_down_H_tilde = zeros(size(x_h_G_down));
x_h_G_down_G_tilde = zeros(size(x_h_G_down));
for i=1:size(x_h_G_down,2)
    x_h_G_down_H_tilde(:,i) = H_tilde(x_h_G_down(:,i)')';
    x_h_G_down_G_tilde(:,i) = G_tilde(x_h_G_down(:,i)')';
end

% vertical downsampling with a factor of 2
x_v_H_down_H = zeros(ceil(size(x_h_H_down,1)/2),size(x_h_H_down,2));
x_v_H_down_G = zeros(ceil(size(x_h_H_down,1)/2),size(x_h_H_down,2));
x_v_G_down_H = zeros(ceil(size(x_h_G_down,1)/2),size(x_h_G_down,2));
x_v_G_down_G = zeros(ceil(size(x_h_G_down,1)/2),size(x_h_G_down,2));
for i=1:size(x_h_H_down,2)
    x_v_H_down_H(:,i) = downsample(x_h_H_down_H_tilde(:,i)',2);
    x_v_H_down_G(:,i) = downsample(x_h_H_down_G_tilde(:,i)',2);
    x_v_G_down_H(:,i) = downsample(x_h_G_down_H_tilde(:,i)',2);
    x_v_G_down_G(:,i) = downsample(x_h_G_down_G_tilde(:,i)',2);
end

% synthesis part
%

% vertical upsampling with a factor of 2
x_v_H_up_H = zeros(2*size(x_v_H_down_H,1),size(x_v_H_down_H,2));
x_v_H_up_G = zeros(2*size(x_v_H_down_G,1),size(x_v_H_down_G,2));
x_v_G_up_H = zeros(2*size(x_v_G_down_H,1),size(x_v_G_down_H,2));
x_v_G_up_G = zeros(2*size(x_v_G_down_G,1),size(x_v_G_down_G,2));
for i=1:size(x_v_H_up_H,2)
    x_v_H_up_H(:,i) = upsample(x_v_H_down_H(:,i),2);
    x_v_H_up_G(:,i) = upsample(x_v_H_down_G(:,i),2);
    x_v_G_up_H(:,i) = upsample(x_v_G_down_H(:,i),2);
    x_v_G_up_G(:,i) = upsample(x_v_G_down_G(:,i),2);
end

% vertical filtering with H(z) and G(z)
x_v_H_H = zeros(size(x_v_H_up_H));
x_v_H_G = zeros(size(x_v_H_up_G));
x_v_G_H = zeros(size(x_v_G_up_H));
x_v_G_G = zeros(size(x_v_G_up_G));
for i=1:size(x_v_H_H,2)
    x_v_H_H(:,i) = H(x_v_H_up_H(:,i)')';
    x_v_H_G(:,i) = H(x_v_H_up_G(:,i)')';
    x_v_G_H(:,i) = G(x_v_G_up_H(:,i)')';
    x_v_G_G(:,i) = G(x_v_G_up_G(:,i)')';
end

comb_H = x_v_H_H + x_v_H_G;
comb_G = x_v_G_H + x_v_G_G;

% horizontal upsampling with a factor of 2
x_h_H_up = zeros(size(comb_H,1),2*size(comb_H,2));
x_h_G_up = zeros(size(comb_G,1),2*size(comb_G,2));
for i=1:size(x_h_H_up,1)
    x_h_H_up(i,:) = upsample(comb_H(i,:),2);
    x_h_G_up(i,:) = upsample(comb_G(i,:),2);
end

% horizontal filtering with H(z) and G(z)
x_h_H = zeros(size(x_h_H_up));
x_h_G = zeros(size(x_h_G_up));
for i=1:size(x_h_H,1)
    x_h_H(i,:) = H(x_h_H_up(i,:));
    x_h_G(i,:) = G(x_h_G_up(i,:));
end




% combine filtered images
x_head = x_h_H + x_h_G;
x_head02 = uint8(x_head);

subplot(1,2,1)
imshow(x_head02)
subplot(1,2,2)
imshow(x)


function x_ext=ext_col_zeros(x,a,b)
    % extend the matrix columns with a zero columns to the left and b zero
    % columns to the right
     x_ext = [ones(size(x,1),a) x ones(size(x,1),b)];
end

function x_ext=ext_row_zeros(x,a,b)
    % extend the matrix rows with a zero rows at the top and b zero
    % rows at the bottom
     x_ext = [zeros(a,size(x,2)); x ;zeros(b,size(x,2))];
end

function x_ext=extention(x,a,b)
    x_ext = ext_col_zeros(x,a,b);
end

function x_conv=my_conv(x,h,h0i)
    % input parameter:
    %   x       signal
    %   h       filter coefficents
    %   h0i     index of filter coeff h0
    x_ext = extention(x,h0i-1,size(h,2)-h0i);
    x_conv = zeros(1,length(x));
    for i=1:length(x)
        for j=1:length(h)
            x_conv(i) = x_conv(i) + h(j)*x_ext(i+j-1);
        end
    end    
end

function y=H(x)
    % lowpass filter H(z) = sqrt(2)/4*(z+2+z^-1)
    h=[sqrt(2)/4,sqrt(2)/2,sqrt(2)/4]; % filter coefficents [h-1,h0,h1]
    y = my_conv(x,h,2);    
end

function y=H_tilde(x)
    % lowpass filter H~(z^-1) = sqrt(2)/8*(-z^-2+2z^-1+6+2z-z^2)
    h=[-sqrt(2)/8,sqrt(2)/4,(3*sqrt(2))/4,sqrt(2)/4,-sqrt(2)/8]; % filter coefficents [h-2,h-1,h0,h1,h2]
    y = my_conv(x,h,3);    
end

function y=G(x)
    % high filter G(z) = sqrt(2)/8*(-z^-3+2z^-2-6z^-1+2-z^1)
    g=[-sqrt(2)/8,sqrt(2)/4,(3*sqrt(2))/4,sqrt(2)/4,-sqrt(2)/8]; % filter coefficents [g-3,g-2,g-1,g0,g1]
    y = my_conv(x,g,4);    
end

function y=G_tilde(x)
    % highpass filter G~(z^-1) = sqrt(2)/4*(-z^-2-2z^-1+1)
    g=[sqrt(2)/4,-sqrt(2)/2,sqrt(2)/4]; % filter coefficents [g-2,g-1,g0]
    y = my_conv(x,g,3);    
end
