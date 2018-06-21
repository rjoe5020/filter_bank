clear all;
x = imread('coins.png');

% analysis part
%
% horizontal filtering with H~(z^-1) and G~(z^-1)
x_h_H_tilde = uint8(zeros(size(x)));
x_h_G_tilde = uint8(zeros(size(x)));
for i=1:size(x,1)
    x_h_H_tilde(i,:) = H_tilde(x(i,:));
    x_h_G_tilde(i,:) = G_tilde(x(i,:));
end

% horizontal downsampling with a factor of 2
x_h_H_down = uint8(zeros(size(x,1),ceil(size(x,2)/2)));
x_h_G_down = uint8(zeros(size(x,1),ceil(size(x,2)/2)));
for i=1:size(x,1)
    x_h_H_down(i,:) = downsample(x_h_H_tilde(i,:),2);
    x_h_G_down(i,:) = downsample(x_h_G_tilde(i,:),2);
end

% vertical filtering with H~(z^-1) and G~(z^-1)
x_v_H_tilde = uint8(zeros(size(x_h_H_down)));
x_v_G_tilde = uint8(zeros(size(x_h_G_down)));
for i=1:size(x_h_H_down,2)
    x_v_H_tilde(:,i) = H_tilde(x_h_H_down(:,i)')';
    x_v_G_tilde(:,i) = G_tilde(x_h_G_down(:,i)')';
end

% vertical downsampling with a factor of 2
x_v_H_down = uint8(zeros(ceil(size(x_h_H_down,1)/2),size(x_h_H_down,2)));
x_v_G_down = uint8(zeros(ceil(size(x_h_G_down,1)/2),size(x_h_G_down,2)));
for i=1:size(x_h_H_down,2)
    x_v_H_down(:,i) = downsample(x_v_H_tilde(:,i)',2);
    x_v_G_down(:,i) = downsample(x_v_G_tilde(:,i)',2);
end

% synthesis part
%
% horizontal upsampling with a factor of 2
if mod(size(x_v_H_down,2), 2) == 0
  x_h_H_up = uint8(zeros(size(x_v_H_down,1),2*size(x_v_H_down,2)));
  x_h_G_up = uint8(zeros(size(x_v_G_down,1),2*size(x_v_G_down,2)));
  for i=1:size(x_v_H_down,1)
      x_h_H_up(i,:) = upsample(x_v_H_down(i,:),2);
      x_h_G_up(i,:) = upsample(x_v_G_down(i,:),2);
  end
else % odd number of columns --> have to remove last column after upsampling to keep correct size
  x_h_H_up = uint8(zeros(size(x_v_H_down,1),2*size(x_v_H_down,2)));
  x_h_G_up = uint8(zeros(size(x_v_G_down,1),2*size(x_v_G_down,2)));
  for i=1:size(x_v_H_down,1)
      x_h_H_up(i,:) = upsample(x_v_H_down(i,:),2);
      x_h_G_up(i,:) = upsample(x_v_G_down(i,:),2);
  end
  x_h_H_up(:,end) = [];
  x_h_G_up(:,end) = [];
end

% horizontal filtering with H(z) and G(z)
x_h_H = uint8(zeros(size(x_h_H_up)));
x_h_G = uint8(zeros(size(x_h_G_up)));
for i=1:size(x_h_H,1)
    x_h_H(i,:) = H(x_h_H_up(i,:));
    x_h_G(i,:) = G(x_h_G_up(i,:));
end

% vertical upsampling with a factor of 2
% if ??
if mod(size(x_h_H,1), 2) == 0
  x_v_H_up = uint8(zeros(2*size(x_h_H,1),size(x_h_H,2)));
  x_v_G_up = uint8(zeros(2*size(x_h_G,1),size(x_h_G,2)));
  for i=1:size(x_h_H,2)
      x_v_H_up(:,i) = upsample(x_h_H(:,i),2);
      x_v_G_up(:,i) = upsample(x_h_G(:,i),2);
  end
else % odd number of columns --> have to remove last column after upsampling to keep correct size
  x_v_H_up = uint8(zeros(2*size(x_h_H,1),size(x_h_H,2)));
  x_v_G_up = uint8(zeros(2*size(x_h_G,1),size(x_h_G,2)));
  for i=1:size(x_h_H,2)
      x_v_H_up(:,i) = upsample(x_h_H(:,i),2);
      x_v_G_up(:,i) = upsample(x_h_G(:,i),2);
  end 
end

% vertical filtering with H(z) and G(z)
x_v_H = uint8(zeros(size(x_v_H_up)));
x_v_G = uint8(zeros(size(x_v_G_up)));
for i=1:size(x_v_H,2)
    x_v_H(:,i) = H(x_v_H_up(:,i)')';
    x_v_G(:,i) = G(x_v_G_up(:,i)')';
end

% combine filtered images
x_head = x_v_H + x_v_G;
subplot(1,2,1)
imshow(x_head)
subplot(1,2,2)
imshow(x)



function x_ext=ext_col_zeros(x,a,b)
    % extend the matrix columns with a zero columns to the left and b zero
    % columns to the right
     x_ext = [zeros(size(x,1),a) x zeros(size(x,1),b)];
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
    %g=[sqrt(2)/8,sqrt(2)/4,-(3*sqrt(2))/4,sqrt(2)/4,sqrt(2)/8];
    y = my_conv(x,g,5);    
end

function y=G_tilde(x)
    % highpass filter G~(z^-1) = sqrt(2)/4*(-z^-2-2z^-1+1)
    g=[1/(2*sqrt(2)),-1/sqrt(2),1/(2*sqrt(2))]; % filter coefficents [g-2,g-1,g0]
    y = my_conv(x,g,2);    
end

