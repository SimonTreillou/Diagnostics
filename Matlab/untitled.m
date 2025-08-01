RGB = imread('/Users/simon/Desktop/test3.png');

[X,MAP] = rgb2ind(RGB,1000);

X=flip(X,2);

colormap summer;
s=pcolor(T3);
set(s,'EdgeColor','none');
colorbar();