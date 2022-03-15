close all;
clear;

w=100e-9;
l=200e-9;
dx=25e-10;
dy=25e-10;
ny=int64(w/dx);
nx=int64(l/dy);

V0=0.1;
cMap=ones(nx,ny);
% Change the conductance for the middle third on each side 
cMap(nx/3:2*nx/3  , 1:ny/3)=1e-15; 
cMap(nx/3:2*nx/3 , 2*ny/3:ny)=1e-15;

% G-matrix, relates the value of a node to all other nodes
G=sparse(nx*ny,nx*ny);
% F-vector, the boundary conditons
F = sparse(nx*ny,1);

% Populate G matrix
% Code from https://github.com/L28E/4700Code/blob/master/CondCode/GetCurrents.m
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            F(n)=V0;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == 1
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nyp = j + 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            rym = (cMap(i, j) + cMap(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end

% Solve for voltage
V = G\F;
for x=1:nx
    for y=1:ny
      n = y + (x - 1) * ny;
      v_surf(x,y)=V(n);
    end      
end 

[X,Y]=meshgrid(linspace(0,l,nx),linspace(0,w,ny));

% Voltage Plot 
figure(1);
surf(X,Y,v_surf','EdgeColor','none');
title('V(x,y)');
ylabel('W');
xlabel('L');

% E field Plot 
[Ex,Ey]=gradient(-v_surf',dx,dy);
figure(2);
quiver(X,Y,Ex,Ey);
title('E(x,y)' );
ylabel('W');
xlabel('L');

