l0 = 1.55e-6;
k0 = 2*pi/l0;

fix = 0;
fiy = 0;

dephase_x = (pi*fix/180);
dephase_y = (pi*fiy/180);

dx = l0/8;
dy = dx;

gamma = 2;

%%%% x & y rango
rX = 2000e-6;
rY = rX;

x = -rX/2:rX/500:rX/2;
y = -rY/2:rY/500:rY/2;

M = length(x);

[X,Y] = meshgrid(x,y);

%%%%%%% rectangular 2D array
dimx = 30e-6;
dimy = 30e-6;

E = zeros(size(X)); % Inicialización correcta

elex = 11; 
eley = 11; % número de elementos

%%%%% configuración del array
%sepx = 1.9*l0;
%sepy = 1.9*l0;
sepx = 148e-6;
sepy = 148e-6;


for ky = -(eley-1)/2 : (eley-1)/2
    for kx = -(elex-1)/2 : (elex-1)/2

        aperture = (abs(X - kx*sepx) < (dimx/2)) .* ...
                   (abs(Y - ky*sepy) < (dimy/2));
               
        E = E + aperture .* exp(1i*(kx*dephase_x + ky*dephase_y));
    end
end

%%%%%% Transformación 2D
FarField = fftshift(fft2(E, gamma*M, gamma*M));

FarField = FarField ./ max(max(abs(FarField)));
FarField = FarField + 1e-9;

fsx = 1/(x(2)-x(1));
sx = linspace(-fsx/2, fsx/2, gamma*M);
angx = 180*real(asin(sx*l0))/pi;

fsy = 1/(y(2)-y(1));
sy = linspace(-fsy/2, fsy/2, gamma*M);
angy = 180*real(asin(sy*l0))/pi;


%%%%%% Gráfica
figure(1)
h = pcolor(angx, angy, abs(FarField).^2);

clim([0 1])
shading flat;
colormap(hot);

xlabel('Angle in plane X');
ylabel('Angle in plane Y');
xlim([-3.5,3.5]);
ylim([-3.5,3.5]);
axis equal;
colorbar

title('Far Field (power)');


%%%%%%% Gráfico 3D %%%%%%%

theta_x = angx * pi / 180;
theta_y = angy * pi / 180;

[THETA_X, THETA_Y] = meshgrid(theta_x, theta_y);

R = abs(FarField).^2;

X_3d = R .* sin(THETA_X);
Y_3d = R .* sin(THETA_Y);

Z_3d = sqrt(max(0, R.^2 - X_3d.^2 - Y_3d.^2));

figure(2)
h3d = mesh(X_3d, Y_3d, Z_3d);

shading interp;       
colormap(jet);        
colorbar;             

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Spherical Far Field Radiation Pattern');

view(-45, 25);        
grid on;              
axis square;          


%%%%%%% GRÁFICA DEL CAMPO CERCANO  %%%%%%%

figure(3)

NearField_Power = abs(E).^2;
NearField_Power = NearField_Power ./ max(NearField_Power(:));

imagesc(x * 1e6, y * 1e6, NearField_Power);

axis xy;                
colormap(parula);        
colorbar;               
axis equal;             

ancho_total_opa = (elex - 1) * sepx + 2 * dimx; 
zoom_limit = (ancho_total_opa / 2) * 1.5 * 1e6; 

xlim([-zoom_limit, zoom_limit]);
ylim([-zoom_limit, zoom_limit]);


title('Near Field (power)');
xlabel('X (\mum)');
ylabel('Y (\mum)');
