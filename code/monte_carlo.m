%=======================================================================
% ELEC 4700 Assignment 3
% Nicholas Ramkhalawansingh

% Part 1
%=======================================================================
clear
close all

m_0=9.10938e-31;        % electron rest mass (kg)
m_n=0.26*m_0;           % electron effective mass (kg)
T=300;                  % Temperature (K)
k_b=1.380649e-23;       % Boltzmann Constant (J/K)
q = 1.60217653e-19;     % Charge of an electron

V_th=sqrt(2*k_b*T/m_n);   % Thermal velocity (m/s)
tau_mn=0.2e-12; % Mean time between collisions 

num_electrons=1000;
num_steps=1000;
num_traces=10; % Increase the number of tracked electrons for MFP accuracy
ymax=100e-9;
xmax=200e-9;
area=xmax*ymax; % Cross sectional area, for the current calculation
dt=4e-15;
P_scat=1-exp(-dt/tau_mn);

% Intoduce an electric field
Ex=0.1/xmax;    % Field is constant throughout material; a function of the voltage and length of the material
Ey=0;
Fx=q*Ex;
Fy=q*Ey;

fprintf("E field: %e V/m\n",sqrt(Ex^2+Ey^2))
fprintf("F on each eletron: %e N\n",sqrt(Fx^2+Fy^2))

% Generate random electron positions
Px=rand(1,num_electrons).*xmax;
Py=rand(1,num_electrons).*ymax;

% Generate random electron velocities (Normal distribution for each component of velocity)
%twiddle=1.25;  
twiddle=1; 
Vx=randn(1,num_electrons)*sqrt(k_b*T/m_n)*twiddle;
Vy=randn(1,num_electrons)*sqrt(k_b*T/m_n)*twiddle;

% Randomly select some electrons to follow
tracked_indices=randperm(num_electrons,num_traces);

% Make vectors to store the paths of those electons, and temperature
X=zeros(num_traces,num_steps);
Y=zeros(num_traces,num_steps);
t=zeros(1,num_steps);
current=zeros(1,num_steps);

% 2D array to track the timesteps where each electron has a collision
collisions=zeros(num_electrons,num_steps);

figure(1)

for k=2:num_steps
    % Update positions
    Px=Px+Vx*dt;
    Py=Py+Vy*dt;
    
    % Scatter electrons
    scat=rand(1,num_electrons)<P_scat;
    Vx(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n)*twiddle;
    Vy(scat)=randn(1,length(Vx(scat)))*sqrt(k_b*T/m_n)*twiddle;    
    
    % Acceleration due to the electric field
    Vx=Vx-Fx/m_n*dt;
    Vy=Vy-Fy/m_n*dt;    
    
    % Electrons leaving lateral bounds come back in to preserve density
    Px(Px<0)=xmax+Px(Px<0);
    Px(Px>xmax)=Px(Px>xmax)-xmax;
    
    % Electrons reflect off upper and lower bounds
    beyond_upper=Py>ymax;
    beyond_lower=Py<0;
    Vy(beyond_lower|beyond_upper)=-Vy(beyond_lower|beyond_upper);
    Py(beyond_lower)=-Py(beyond_lower);
    Py(beyond_upper)=-Py(beyond_upper)+2*ymax;    
    
    % Update the tracked electrons and temperature 
    t(k)=t(k-1)+dt;
    X(:,k)=Px(tracked_indices);
    Y(:,k)=Py(tracked_indices);    
    current(k)=mean(Vx) *q *10e19 ... % Current Density 
               *area; % Cross sectional area
    
    % Record the time steps where the electrons scatter
    %collisions(scat|beyond_upper|beyond_lower,k)=1;
    collisions(scat,k)=1;
    
    % Plot electron trajectories        
    plot(X(1,1:k),Y(1,1:k),".",X(2,1:k),Y(2,1:k),".",X(3,1:k),Y(3,1:k),".",...
        X(4,1:k),Y(4,1:k),".",X(5,1:k),Y(5,1:k),".",'MarkerSize',4) % hide the discontinuities by using dots :)
       
   
    pause(0.00001)    
end
title("1C: Electron Trajectories")
ylabel("Y (m)")
xlabel("X (m)")

% Plot current
figure(2)
plot(t,abs(current))
title("1D: Current vs. Time")
ylabel("Current (A)")
xlabel("Time (s)")

% Plot temperature
%https://www.mathworks.com/help/matlab/ref/griddata.html
figure(3);
subplot(1,2,1);
temp=(Vx(:).^2+Vy(:).^2)*m_n/k_b/2;
[xq,yq] = meshgrid(0:0.05*xmax:xmax,0:0.05*ymax:ymax);
vq = griddata(Px,Py,temp',xq,yq);
s=mesh(xq,yq,vq);
s.FaceColor='interp';
c=colorbar;
c.Label.String="Temperature (k)";
title("Temperature Map")
ylabel("Y (m)")
xlabel("X (m)")

% Plot density map
subplot(1,2,2);
binscatter(Px,Py,20)
title("Electron Density Map")
ylabel("Y (m)")
xlabel("X (m)")
axis([0 xmax 0 ymax])
colormap('parula')
