close all; clearvars; clc;

load('synthdata.mat')

m = Depth(:);
Image = Image; %#ok
fd = Dparam.fp;
fc = Cparam.fp;
cx = Cparam.cx;
cy = Cparam.cy;
R = R; %#oksiz
T = T; %#ok

%% Task 1

size_pic_d = size(Depth);
size_pic = size(Image);
[xxi, yyi] = meshgrid(1:size_pic_d(2), 1:size_pic_d(1));

%Fixing to the optical origen
xxi = xxi - size_pic_d(2)/2;
yyi = yyi - size_pic_d(1)/2;

xx = xxi(:);
yy = yyi(:);

zz = fd*ones(size(xx));

z_i = fd*(m./(sqrt(xx.^2 + yy.^2 + zz.^2)));

up = xx.*z_i;
vp = yy.*z_i;

M = [fd,  0,  0;
     0,  fd,  0;
     0,  0,  1];

aux = M\[up'; vp'; z_i'];

x = aux(1, :);
y = aux(2, :);
z = aux(3, :);

figure(1)
scatter3(x, y, z, 2, z)
zlim([0, max(z)])
set(gca, 'Zdir', 'reverse')
set(gca,'YDir','reverse')
view([0 90])
title('Task 1: Point cloud in global (x,y,z) space')

hold on
plot(0,0,'rx', 'MarkerSize', 20, 'LineWidth', 2)
hold off

%% Task 2

M2 = [fc,  0,  cx;
      0,  fc,  cy;
      00,  0,  1];

aux2ax = [inv(R), -T]*[x; y; z; ones(size(x))];
aux2 = M2*aux2ax;

u = aux2(1, :)./aux2(3, :);
v = aux2(2, :)./aux2(3, :);
z2 = aux2(3, :);

figure(2)
scatter(u, v, 2, z2)
axis([0 2*cx 0 2*cy])
set(gca,'YDir','reverse')
title('Task 2: Global depth points projected on image plane of the color camera (IRREGULAR)')

%% Task 3

F = scatteredInterpolant(double(u'), double(v'), double(z2'), 'nearest', 'nearest');

[xq, yq] = meshgrid(1:2*cx, 1:2*cy);

vq = F(xq, yq);

%% Task 4

figure(3)
surf(xq, yq, vq, Image);
axis([0 2*cx 0 2*cy])
set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
view([40 20])
shading interp
title('Task 3+4: Global depth points projected on image plane of the color camera (RESAMPLED)')

%% Task 5

figure(4)
h = surf(xq, yq, vq, Image);
remart3D(h, size_pic_d);

set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
view([50 20])
shading interp
title('Task 5: Removing 3D model edge artifacts')

%% Task 6

vr_6 = reshape(z', size_pic_d(1:2));
u_6 = reshape(u', size_pic_d(1:2));
v_6 = reshape(v', size_pic_d(1:2));

vq_6 = interp2(vr_6, xq, yq);

Col = zeros(size(Image));
for o = 1:3 
    Col(:, :, o) = interp2(Image(:, :, o), u_6, v_6);
end

figure(5)
h2 = surf(vq_6, Col);
remart3D(h2, size_pic_d);

set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
view([40 30])
shading interp

title('Task 6: Global depth points projected on image plane of the color camera (INTERPOLATION)')

%% Task 7
zau = z2(:);
au = z2(:);

IDX = knnsearch([u(:), v(:)], [u(:), v(:)], 'K', 3);
new = {0};

for o = 1:length(u(:))
    ner = IDX(o, :);
    val = min(au(ner));
    
    new{o} = ner(au(ner) > 1.1*val);
    
    zau(ner) = val;
end

F = scatteredInterpolant(double(u'), double(v'), double(zau), 'nearest', 'nearest');
vq_7 = F(xq, yq);

figure(6)
h3 = surf(xq, yq, vq_7, Image);
remart3D(h3, size_pic_d);
set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
view([50 20])
shading interp
title('Task 7: Removing 3D model edge artifacts (Z-Buffering)')

%% Task 8

vocl = horzcat(new{:})';

xq(vocl) = NaN;
yq(vocl) = NaN;
u_6(vocl) = NaN;
v_6(vocl) = NaN;

vq_8 = interp2(vr_6, xq, yq);

Col = zeros(size(Image));
for o = 1:3
    Col(:, :, o) = interp2(Image(:, :, o), u_6, v_6);
end

figure(7)
h4 = surf(vq_8, Col);
remart3D(h4, size_pic_d);

set(gca,'ZDir','reverse')
set(gca,'YDir','reverse')
view([50 20])
shading interp
title('Task 8: Removing 3D model edge artifacts (Z-Buffering / INTERPOLATION)')