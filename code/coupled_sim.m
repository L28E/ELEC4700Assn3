function [mean_current]=coupled_sim(question,V0,channel_width)
% Part 2c,3a,3b
% Usage: 
% coupled_sim("2C",0.1,0.2e-7)
% coupled_sim("3A",0.8,0.2e-7)
% coupled_sim("3B",0.1,0.1e-7)
       
    type=0;                 % 0 for specular, anything else for diffusive 

    m_0=9.10938e-31;        % electron rest mass (kg)
    m_n=0.26*m_0;           % electron effective mass (kg)
    T=300;                  % Temperature (K)
    k_b=1.380649e-23;       % Boltzmann Constant (J/K)
    tau_mn=0.2e-12;         % Mean time between collisions 
    q = 1.60217653e-19;     % Charge of an electron

    num_electrons=1000;
    num_steps=1000;
    num_traces=5;
    ymax=100e-9;
    xmax=200e-9;    
    dt=4e-15;
    P_scat=1-exp(-dt/tau_mn);

    % Generate the boxes
    box_x=0.8e-7;
    box_width=0.4e-7;
    box_right=box_x+box_width;
    
    box_height=(ymax-channel_width)/2;
    box_y=0;
    box_top=box_y+box_height;
    box_bottom=ymax-box_height;

    % Generate random electron positions
    Px=rand(1,num_electrons).*xmax;
    Py=rand(1,num_electrons).*ymax;

    % Remove electrons from the boxes
    oob=get_oob(Px,Py,box_x,box_right,box_top,box_bottom);
    [M,I]=get_nearest_bound(num_electrons,Px,Py,oob,box_x,box_right,box_top,box_bottom);
    [Px,Py]=remove_oob(Px,Py,oob,box_x,box_right,box_top,box_bottom,M,I);

    % Generate random electron velocities (Normal distribution for each component of velocity)
    Vx=randn(1,num_electrons)*sqrt(k_b*T/m_n);
    Vy=randn(1,num_electrons)*sqrt(k_b*T/m_n);

    % Randomly select some electrons to follow
    tracked_indices=randperm(num_electrons,num_traces);

    % Make vectors to store the paths of those electons, and temperature
    X=zeros(num_traces,num_steps);
    Y=zeros(num_traces,num_steps);
    t=zeros(1,num_steps);
    current=zeros(1,num_steps);

    % 2D array to track the timesteps where each electron has a collision
    collisions=zeros(num_electrons,num_steps);

    % Get the Electric field using finite difference method
    dx=25e-10;
    dy=dx;
    [Ex,Ey,nx,ny]=bottleneck_E_field(false,V0,ymax,xmax,dx,dy,box_x,box_right,box_top,box_bottom);

    if question=="2C"
        figure(1);
    end
    
    for k=2:num_steps
        % Update positions
        Px=Px+Vx*dt;
        Py=Py+Vy*dt;

        % Scatter electrons
        scat=rand(1,num_electrons)<P_scat;
        Vx(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n);
        Vy(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n);    

        % TODO: Try to optimize this    
        %=================================================    
        % Force is a function of the E-field, which is a function of space because of the bottleneck.
        % For all electrons in a dx*dy E-field "bin", apply the force of that "bin"       

        for x=1:nx
           for y=1:ny
               Vx(Px>(double(x)-1)*dx & Px<double(x)*dx & Py>(double(y)-1)*dy & Py<double(y)*dy)=Vx(Px>(double(x)-1)*dx & Px<double(x)*dx & Py>(double(y)-1)*dy & Py<double(y)*dy)-q*Ex(y,x)/m_n*dt;
               Vy(Px>(double(x)-1)*dx & Px<double(x)*dx & Py>(double(y)-1)*dy & Py<double(y)*dy)=Vy(Px>(double(x)-1)*dx & Px<double(x)*dx & Py>(double(y)-1)*dy & Py<double(y)*dy)-q*Ey(y,x)/m_n*dt;            
           end      
        end    
        %=================================================

        % Electrons leaving lateral bounds come back in to preserve density
        Px(Px<0)=xmax+Px(Px<0);
        Px(Px>xmax)=Px(Px>xmax)-xmax;

        % Electrons reflect off upper and lower bounds
        beyond_upper=Py>ymax;
        beyond_lower=Py<0;
        Vy(beyond_lower|beyond_upper)=-Vy(beyond_lower|beyond_upper);
        Py(beyond_lower)=-Py(beyond_lower);
        Py(beyond_upper)=-Py(beyond_upper)+2*ymax;  

        % Logically index electrons which are out of bounds
        oob=get_oob(Px,Py,box_x,box_right,box_top,box_bottom);    

        % Determine which bound each oob electron is closest to
        [M,I]=get_nearest_bound(num_electrons,Px,Py,oob,box_x,box_right,box_top,box_bottom);  

        hzn_flip=oob & M'~=0 & (I'==1|I'==2); % For oob electrons nearest the left or right edges, flip their horizontal velocity
        vert_flip=oob & M'~=0 & (I'==3|I'==4); % For oob electrons nearest the top boundary, flip their vertical velocity

        % Reflect & remove oob electrons from boxes
        if type==0          
            Vx(hzn_flip)=-Vx(hzn_flip);          
            Vy(vert_flip)=-Vy(vert_flip);
        else
            Vx(hzn_flip)=-sign(Vx(hzn_flip)).*abs(randn(1,nnz(hzn_flip))*sqrt(k_b*T/m_n));        
            Vy(vert_flip)=-sign(Vy(vert_flip)).*abs(randn(1,nnz(vert_flip))*sqrt(k_b*T/m_n));
        end
        [Px,Py]=remove_oob(Px,Py,oob,box_x,box_right,box_top,box_bottom,M,I);   

        % Update the tracked electrons and temperature 
        t(k)=t(k-1)+dt;
        X(:,k)=Px(tracked_indices);
        Y(:,k)=Py(tracked_indices);    
        current(k)=mean(Vx) *q *10e19 ... % Current Density 
                   *ymax; % Cross sectional area

        % Record the time steps where the electrons scatter
        %collisions(scat|beyond_upper|beyond_lower,k)=1;
        collisions(scat,k)=1;

        fprintf('%2.1f%%\n',k/num_steps*100)
        if question=="2C"        
            % Plot electron trajectories   
            plot(X(1,1:k),Y(1,1:k),".",X(2,1:k),Y(2,1:k),".",X(3,1:k),Y(3,1:k),".",...
                X(4,1:k),Y(4,1:k),".",X(5,1:k),Y(5,1:k),".",...
                'MarkerSize',4) % hide the discontinuities by using dots :)
            hold on;

            rectangle('Position',[box_x box_y box_width box_height]);
            rectangle('Position',[box_x box_bottom box_width box_height]);
            hold off;
            
            pause(0.00001)        
        end
    end
    
    if question=="2C"
        title("2C: Electron Trajectories")
        ylabel("Y (m)")
        xlabel("X (m)")
    end
    
    % Plot density map
    if question=="3A"
        figure();
        binscatter(Px,Py,24)
        title("3a: Electron Density Map")
        ylabel("Y (m)")
        xlabel("X (m)")
        axis([0 xmax 0 ymax])
        colormap('parula')
    end
    
    if question=="3B"
        % Plot current
        figure();
        plot(t,abs(current))
        title("Current vs. Time")
        ylabel("Current (A)")
        xlabel("Time (s)")
        mean_current=mean(abs(current));
    end
    
end
