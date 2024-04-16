%% Gradient acoustic resonance visualiser
% Reads in simulator gradient files and outputs waveforms & FFT
% PJL 25/01/22
clear all
startdir=pwd;

%% Folder for simulation files
sim_name = "./sims/cones_test";
cd(sim_name);

%% Simulation settings
tw_len_ms   = 100;                 % duration of sliding time window (ms)  
update_rate_ms = tw_len_ms;        % how far to step the sliding time window per frame (ms)

%% Scanner settings - *GET THESE FROM IMPRINTD FOR YOUR SYSTEM*
TS = 1e-5;             % Gradient raster rate (s)
SMAX = 200;            % Slew rate limit (T/m/s)
GMAX = 80;             % Max grad amplitude (mT/m)
forbidden_frq=[550,1100];       % Forbidden frequencies (TERRA)
forbidden_bw=[100,300];         % Forbidden frequency ranges (TERRA)
min_frq=forbidden_frq-forbidden_bw./2;
max_frq=forbidden_frq+forbidden_bw./2;


%% Read in X Gradient

vf_txt = split(regexp(fileread(dir("./*_GRX.dsv").name),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
temp_mat = readmatrix(dir("./*_GRX.dsv").name, "FileType", "text");

XG=[];
XG(1)=0;
midx=1;
gidx=1;

while midx<(size(temp_mat,1))
    midx=midx+1;
    gidx=gidx+1;
    XG(gidx)=temp_mat(midx);
    if temp_mat(midx)==temp_mat(midx-1)
        XG(gidx+1:gidx+temp_mat(midx+1)+1)=XG(gidx);
        gidx=gidx+temp_mat(midx+1)+1;
        XG(gidx)=temp_mat(midx+2);
        midx=midx+2;
    end
end

XG=cumsum(XG(1:end-1))./vert_scl;

%% Read in Y Gradient

vf_txt = split(regexp(fileread(dir("./*_GRY.dsv").name),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
temp_mat = readmatrix(dir("./*_GRY.dsv").name, "FileType", "text");

YG=[];
YG(1)=0;
midx=1;
gidx=1;

while midx<(size(temp_mat,1)-1)
    midx=midx+1;
    gidx=gidx+1;
    YG(gidx)=temp_mat(midx);
    if temp_mat(midx)==temp_mat(midx-1)
        YG(gidx+1:gidx+temp_mat(midx+1)+1)=YG(gidx);
        gidx=gidx+temp_mat(midx+1)+1;
        YG(gidx)=temp_mat(midx+2);
        midx=midx+2;
    end
end

YG=cumsum(YG(1:end-1))./vert_scl;

%% Read in Z Gradient

vf_txt = split(regexp(fileread(dir("./*_GRZ.dsv").name),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
temp_mat = readmatrix(dir("./*_GRZ.dsv").name, "FileType", "text");

ZG=[];
ZG(1)=0;
midx=1;
gidx=1;

while midx<(size(temp_mat,1)-1)
    midx=midx+1;
    gidx=gidx+1;
    ZG(gidx)=temp_mat(midx);
    if temp_mat(midx)==temp_mat(midx-1)
        ZG(gidx+1:gidx+temp_mat(midx+1)+1)=ZG(gidx);
        gidx=gidx+temp_mat(midx+1)+1;
        ZG(gidx)=temp_mat(midx+2);
        midx=midx+2;
    end
end

ZG=cumsum(ZG(1:end-1))./vert_scl;



%% Visualise the sliding time window and FFT
cd(startdir)
tw_len_us= round(tw_len_ms/1000/TS);     % correct for gradient raster rate
update_rate= round(update_rate_ms/1000/TS);   % correct for gradient raster rate  


for tw_no=1:update_rate:length(XG)-tw_len_us

    XG_int=XG(tw_no:(tw_no+tw_len_us-1));
    YG_int=YG(tw_no:(tw_no+tw_len_us-1));
    ZG_int=ZG(tw_no:(tw_no+tw_len_us-1));
    
    h=figure(1);
    
    figure(1) 
    subplot(2,1,1)

    plot((0:tw_len_us-1)*TS*1000,XG_int,"k")
    hold on
    plot((0:tw_len_us-1)*TS*1000,YG_int,"b")
    plot((0:tw_len_us-1)*TS*1000,ZG_int,"r")
    plot([0,(tw_len_us-1)*TS*1000],[GMAX,GMAX],":k")
    plot([0,(tw_len_us-1)*TS*1000],[-GMAX,-GMAX],":k")
    hold off
    xlim([0 tw_len_us*TS*1000])
    ylim([-GMAX*1.2 GMAX*1.2])
    xlabel("Time (ms)","Interpreter","latex","FontSize",12)
    ylabel("Gradient amp. (mT/m)","Interpreter","latex","FontSize",12)

    set(gca,"TickLabelInterpreter","latex","FontSize",12)
    title(strcat("X-Y, Z waveforms during",{' '},num2str(tw_len_us*TS*1e3),"ms sliding window"),"Interpreter","latex","FontSize",14)
    
    
    subplot(2,1,2)
    stem(min_frq(1):max_frq(1),repmat(30,1,max_frq(1)-min_frq(1)+1),":","Color",[0.5,0.5,0.5],"LineWidth",0.2)
    hold on
    stem(min_frq(2):max_frq(2),repmat(30,1,max_frq(2)-min_frq(2)+1),":","Color",[0.5,0.5,0.5],"LineWidth",0.2)

    plot((0:(tw_len_us-1))/TS/tw_len_us,abs(fft(XG_int)./tw_len_us),"-k")
    plot((0:(tw_len_us-1))/TS/tw_len_us,abs(fft(YG_int)./tw_len_us),"-b")
    plot((0:(tw_len_us-1))/TS/tw_len_us,abs(fft(ZG_int)./tw_len_us),"-r")
    xlim([0 3000])
    ylim([0 20])
    hold off
    xlabel("Freq. (Hz)","Interpreter","latex")
    set(gcf,"color","w")
    set(gca,"TickLabelInterpreter","latex","FontSize",12)
    title(strcat("FFT across",{' '},num2str(tw_len_us*TS*1e3),"ms window"),"Interpreter","latex","FontSize",14)
    

    drawnow;
    

end
