function [shots, data] = HLT11_shots()

    %% w46
    ii=1;
    day1 = [64897]; %modes not so clear
    data(ii).shot = day1(1); data(ii).t = [0., 0.]; ii=ii+1;
    day2 = [64935]; % modes not so clear
    data(ii).shot = day2(1); data(ii).t = [0., 0.]; ii=ii+1;
    day3 = [64946, 64949];
    data(ii).shot = day3(1); data(ii).t = [0., 0.]; ii=ii+1;
    data(ii).shot = day3(2); data(ii).t = [0.905, 1.1]; ii=ii+1;
    w46 = [day1, day2, day3];


    %% w47
    day1 = [64994, 64996, 65003, 65021];
    data(ii).shot = day1(1); data(ii).t = [0.702, 0.87];ii=ii+1;
    data(ii).shot = day1(2); data(ii).t = [0.73, 1.];ii=ii+1;
    data(ii).shot = day1(3); data(ii).t = [0., 0.]; ii=ii+1;%modes are strange for 65003
    data(ii).shot = day1(4); data(ii).t = [0., 0.]; ii=ii+1;%modes are strange for 65021

    day2 = [65040, 65050, 65051, 65052, 65055, 65057];
    data(ii).shot = day2(1); data(ii).t = [0.775, 0.92];ii=ii+1;
    data(ii).shot = day2(2); data(ii).t = [0.775, 0.92];ii=ii+1;
    data(ii).shot = day2(3); data(ii).t = [0.5, 0.8];ii=ii+1; %weak modes
    data(ii).shot = day2(4); data(ii).t = [0.55, 0.7];ii=ii+1; %only pos. triangularity case
    data(ii).shot = day2(5); data(ii).t = [0.55, 0.7];ii=ii+1; 
    data(ii).shot = day2(6); data(ii).t = [0.55, 0.7];ii=ii+1; %only pos. triangularity case
    day3 = [65089, 65092, 65094];
    data(ii).shot = day3(1); data(ii).t = [0.45, 1.1];ii=ii+1;
    data(ii).shot = day3(2); data(ii).t = [0.8, 1.1];ii=ii+1;
    data(ii).shot = day3(3); data(ii).t = [0.55, 1.2];ii=ii+1;

    w47 = [day1, day2, day3];

    shots=[w46, w47];
end