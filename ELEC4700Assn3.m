addpath('./code')

monte_carlo
disp("Done Part 1. Press any key to continue...");
pause;
clear;
close all;

bottleneck_E_field(true,0.1,100e-9,200e-9,25e-10,25e-10,0.8e-7,1.2e-07,4e-8,6e-8);
disp("Done Part 2 a,b. Press any key to continue...");
pause;
close all;

coupled_sim("2C",0.1,0.2e-7);
disp("Done Part 2 c. Press any key to continue...");
pause;
close all;

coupled_sim("3A",0.8,0.2e-7);
disp("Done Part 3 a. Press any key to continue...");
disp("WARNING: The next one takes a long time");
pause;
close all;


w=0.1e-7:0.05e-7:0.75e-7;
avg_currents=zeros(1,length(w));
for k=1:length(w)
    avg_currents(k)=coupled_sim("3B",0.1,w(k));
end

plot(avg_currents,w);
title('3b: Average Current vs. Channel Width');
ylabel('Average Current (A)');
xlabel('Channel Width (m)');   

disp("Done Part 3 b. Press any key to continue...");
pause;
clear;
close all;
