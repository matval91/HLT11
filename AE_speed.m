function [va, freq]=AE_speed(shot, rhopos, liuqe,trange, varargin)

% function to plot the spectrogram, AE speed and other stuff
% shot = shot number
% rhopos = position where to compute the AE frequency
% varargin can be one of the following: z, delta, polangle, ne

% example
% AE_speed(shot, 0.9, 'z','delta', 'polangle2')


set(0,'DefaultlineLineWidth',2); % do use width of 2 by default set in set_defaults_matlab
set(0,'defaultAxesFontSize',16);
colors = ['w', 'y', 'm', 'g'];

%% input arguments
if nargin<2; error('Shot and rhopos are mandatory'); end
if ~exist('liuqe','var') || isempty(liuqe); liuqe=0;end
if ~exist('trange','var') || isempty(trange); trange=[0.1, 2];end

opt_args = ["z", "delta", "polangle2", "nel", "ecpower"]; 
offset = 300; %to plot it good in the spec
nn = [500, 200, 2, 2e-18, 1./4]; %to plot it good in the spectrogram
ii=0;
for i=1:length(opt_args)
    if find(contains(varargin,opt_args(i)))
        disp(opt_args(i))
        ii=ii+1;
        if opt_args(i)=='polangle2'
            data_toread{ii} = '\tcv_shot::top.results.toray.input:theta_launch';
        elseif opt_args(i)=='z'
            data_toread{ii} = 'z_axis';
        elseif opt_args(i) == 'ecpower'
            data_toread{ii} = 'powers';
        else
            data_toread{ii} = opt_args(i);
        end
    else
        nn(i) = 0;
    end
end
normalization = nn(nn~=0);
figure();
spec_figure = read_spectrogram(shot, trange);
colormap hot; title(num2str(shot))
axis xy; caxis([0.5 2]); xlabel('t (s)'); ylabel('f (kHz)');
ylim([30, 400]); xlim(trange)
hold on; 
for i=1:length(rhopos)
    disp(rhopos(i))
    [va, freq] = alfven_speed(shot, rhopos(i), liuqe, trange);
    plot(freq.t, freq.data*1e-3, 'c');
    txt = sprintf('rho_p=%.2f', rhopos(i));
    text(max(freq.t)*0.8, mean(freq.data*1e-3)-80, txt, 'Color', 'c', 'FontSize', 18);
    txt = sprintf('q~%.2f', freq.q);
    text(max(freq.t)*0.8, mean(freq.data*1e-3)-100, txt, 'Color', 'c', 'FontSize', 18);
end
if exist('data_toread') && ~isempty(data_toread)
    for i=1:length(data_toread)
        dd=gdat_tcv(shot, data_toread{i},'time_out',freq.t, 'liuqe',liuqe);
        if data_toread{i}=="\tcv_shot::top.results.toray.input:theta_launch"
            dd.data=dd.data(2,:);
        elseif data_toread{i} == "powers"
            dd=gdat_tcv(shot, data_toread{i}, 'liuqe',liuqe);
            dd.data=dd.ec.data(:,10)*0.0001 -offset;
            dd.t = dd.ec.t;
        end
        data = dd.data*normalization(i)+offset;
        plot(dd.t, data,[colors(i),'-']);
        
        if data_toread{i}=="z_axis" || data_toread{i}=="delta"
            [~, ind] = min(abs(max(dd.data)-dd.data)>0.);
            txt = sprintf('%s=%.2f', data_toread{i}, dd.data(ind));
            text(dd.t(ind), data(ind)+20,txt,  'Color', colors(i), 'FontSize', 12)
        
            [~, ind] = min(abs(dd.data-min(dd.data)>0.));
            txt = sprintf('%s=%.2f', data_toread{i}, dd.data(ind));
            text(dd.t(ind), data(ind)-20,txt,  'Color', colors(i), 'FontSize', 12)
        elseif data_toread{i}=="\tcv_shot::top.results.toray.input:theta_launch"
            [~, ind] = min(abs(max(dd.data)-dd.data)>0.);
            txt = sprintf('L2 pol=%.2f', dd.data(ind));
            text(dd.t(ind), data(ind)+20,txt,  'Color', colors(i), 'FontSize', 12)
        
            [~, ind] = min(abs(dd.data-min(dd.data)>0.));
            txt = sprintf('L2 pol=%.2f', dd.data(ind));
            text(dd.t(ind), data(ind)-20,txt,  'Color', colors(i), 'FontSize', 12)
        else
            txt = sprintf('nel=%.2e', mean(dd.data));
            text(dd.t(end)-0.1, mean(data)+50,txt,  'Color', colors(i), 'FontSize', 12)            
        end
    end
else
    disp('No data to plot')
end
hold off;
end

function [va, freq] = alfven_speed(shot, rhopos, liuqe, trange)
    %% alfven frequency speed
    nes=gdat(shot, 'ne_rho'); ne=nes.fit.data;
    Zeff = gdat(shot, 'zeff'); Zeff=Zeff.data;
    BT = gdat_tcv(shot, 'b0', 'time_out', nes.fit.t); BT=BT.data;
    Zc=6; Zi=1;
    Mi=2*1.6e-27; Mc=12*1.6e-27;
    ni     = ne'.*(Zc-Zeff)/(Zc-Zi)/Zi;
    nc     = ne'.*(Zeff-Zi)/(Zc-Zi)/Zc;
    Va     = abs(BT)./sqrt(4*pi*1e-7*(Mi*ni+Mc*nc)); %time, rho

    qs=gdat_tcv(shot, 'q_rho','time_out',nes.fit.t, 'liuqe', liuqe); q=qs.data;
    R0=gdat_tcv(shot, 'r_axis','time_out',nes.fit.t,  'liuqe', liuqe);R0=R0.data;

    ind_Va=iround(nes.fit.x, rhopos);
    ind_q = iround(qs.x, rhopos);
    frequence=Va(:,ind_Va)./(4*pi*q(ind_q,:)'.*R0);
    iind=(nes.fit.t>min(trange) & nes.fit.t<max(trange));
    freq.t = nes.fit.t(iind); freq.data = frequence(iind);freq.q = mean(q(ind_q,:));
    va.t = nes.fit.t(iind); va.data = Va(iind);
end

function spec_figure = read_spectrogram(shot, trange)
%% LTCC magnetic coils

%% colormaps
load /macii/karpusho/Work/2019_MST1_HLT11/AK_ColMaps.mat
%colormap('hot');
colormap(AK_ColMapCold2Hot)
% max fluctuation amplitude: 1(pol), 2 (rad), 4(pol), 5(rad)
% % mni fluctuation amplitude: 7(pol), 8(rad)
% signal = '\ATLAS::DT132_LTCC_001:CHANNEL_002'; %OK
% %signal = '\ATLAS::DT132_LTCC_001:CHANNEL_005';
% %signal = '\ATLAS::DT132_LTCC_001:CHANNEL_008';
% %signal = '\ATLAS::DT132_LTCC_001:CHANNEL_002';
% %signal = '\ATLAS::DT132_LTCC_001:CHANNEL_002';
% 
% if max(trange)-min(trange)<1.
%     n=1;
% else
%     n=4;
% end
% nfft = n*1024;
% 
% 
% % read the data:
% disp(['reading #' num2str(shot) ' ' char(signal) ', please wait ...']);
% tic % start watch to see how much time this routine takes
% %  trace=SF2ML(signal,diag,shot);
%   server = 'tcvdata.epfl.ch';
%   mdsconnect(server);
%   mdsopen('tcv_shot',shot);
%   mag=mdsvalue(signal);
%   tmag=mdsvalue(['dim_of(' signal ')']);
%   mdsclose();
% toc
% if length(mag) > 0
%     ndat = length(tmag);
%     t0 = tmag(1);
%     dt = (tmag(ndat)-t0)/ndat;
%     fsamp = 1 / dt;
% else
%     error('cannot read LPF data');
% end
% 
% % Start of the display section:
% [spec,freq_spec,time_spec] = specgram(mag(:,1),nfft,fsamp);
% time_spec = time_spec + t0;
% freq_spec = freq_spec/1.e3;
% 
% %spec      = sqrt(abs(spec));
% spec      = log10(abs(spec));
% 
% 
% sigma=0.9;
% spec_smooth = imgaussfilt(spec,sigma);
% 
% figure();
% 
% imagesc(time_spec,freq_spec,spec_smooth);
% xl=trange;
% if max(trange)>max(time_spec)
%     xl=[min(time_spec), max(time_spec)];end
% xlim(xl),xlabel('Time (s)');
% ylim([10,400]),ylabel('Frequency (kHz)'), title(['TCV #' num2str(shot) ' ' char(signal)])
% axis xy;                   % make axis go from 0 to max
% set(gca,'fontsize',10); % increase font size
% set(gca,'FontName','arial'); % chnage font type
% h = colorbar;              % add the colorbar
% %set(h,'fontsize',8);
% set(h,'FontName','arial'); % chnage font type
% colormap hot;
% ma = max(caxis)-1;
% mi = min(caxis)+1;
% mi=-0.2; ma=2;
% caxis([mi ma]); % reset plot range
% 
% spec_figure = gca;
% 
% close all;

D = LTCC_spec(shot);
close all;

figure(); 
D.T = D.T(2); D.T = D.T{1};
D.selT = D.selT(2); D.selT = D.selT{1};
D.F = D.F(2); D.F = D.F{1};
D.selF = D.selF(2); D.selF = D.selF{1};
D.Slog=D.Slog(2); D.Slog = D.Slog{1};
imagesc(D.T(D.selT),D.F(D.selF,1),D.Slog(D.selF,D.selT));

spec_figure = gca;

end