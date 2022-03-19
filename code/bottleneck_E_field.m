function [Ex,Ey,nx,ny] = bottleneck_E_field(main,V0,w,l,dx,dy,bottle_left,bottle_right,bottle_top,bottle_bottom)
% Part 2a,b
% Computes the E field of the bottleneck scenario.
% Usage: bottleneck_E_field(true,0.1,100e-9,200e-9,25e-10,25e-10,0.8e-7,1.2e-07,4e-8,6e-8) 

    ny=int64(w/dx);
    nx=int64(l/dy);
    
    % Change the conductance for the bottleneck on each side 
    cMap=ones(nx,ny);
    cMap(int64(bottle_left/dx):int64(bottle_right/dx) , 1:int64(bottle_top/dy))=1e-2; 
    cMap(int64(bottle_left/dx):int64(bottle_right/dx) , int64(bottle_bottom/dy):ny)=1e-2;

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
    v_surf=zeros(nx,ny);
    V = G\F;    
    for x=1:nx
        for y=1:ny
          n = y + (x - 1) * ny;
          v_surf(x,y)=V(n);
        end      
    end
    
    [Ex,Ey]=gradient(-v_surf',dx,dy);    
    
    if main==true
        [X,Y]=meshgrid(linspace(0,l,nx),linspace(0,w,ny));

        % Voltage Plot 
        figure(1);
        surf(X,Y,v_surf','EdgeColor','none');
        title('2a: V(x,y)');
        ylabel('W (m)');
        xlabel('L (m)');
        zlabel('Volts (V)')
        c=colorbar;
        c.Label.String="Volts (V)";

        % E field Plot
        figure(2);
        quiver(X,Y,Ex,Ey);
        title('2b: E(x,y)' );
        ylabel('W (m)');
        xlabel('L (m)');        
    end    
   
end

