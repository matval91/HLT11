function summary_plot_wech(shotin, save)
% For a given shot(s), plots some useful traces
%
% CALL:     summary_plot_os(shot,tm,tM)
%           summary_plot_os(shot);
%           summary_plot_os(shot,[],[],0); % do not print (default)
%
% doprint=0 by default in _os case
% doprint = 0: do not print any
%         = 1: make png and save to Terra as original summary_plot
%         = 2: make png to local folder
%         = -1, -2; as 1, 2 but close figure after each phase (good with multiple shots)
%
% freq_max  max freq [kHz] for spectrogram (default 20)
%
% Benoit Labit - SPC/EPFL - Nov. 2017

ldisplay = matlab.ui.internal.isFigureShowEnabled;
set(0,'Units','Pixels');
sppi  = get(0,'ScreenPixelsPerInch');

if nargin<1; error('Shot number is mandatory');end

%%
set_defaults_matlab % to have color order from colos irrespective of matlab version
set(0,'DefaultlineLineWidth',2); % do use width of 2 by default set in set_defaults_matlab
set(0,'defaultAxesFontSize',16);
%% loop for reading shotsset
for ishot=1:length(shotin)
  %% read ip, ne, z, triang, ec angles
  ip=gdat(shotin(ishot), 'ip'); ne = gdat(shotin(ishot), 'nel'); 
  z = gdat(shotin(ishot), 'z_axis'); delta = gdat(shotin(ishot), 'delta');
  powers = gdat(shotin(ishot), 'powers');
%   b0 = gdat(shotin(ishot), 'b0');
  
  theta_angle_launcher=gdat(shotin(ishot), '\tcv_shot::top.results.toray.input:theta_launch', 'doplot',0);
  phi_angle_launcher=gdat(shotin(ishot), '\tcv_shot::top.results.toray.input:phi_launch', 'doplot',0);
  
  %% figure with (ip,ne), (powers), (angles_theta), (z,delta)
  figure('Position', [10, 10, 1000, 800])
  ax1=subplot(2,2,1); hold on; title(num2str(shotin(ishot)));
  plot(ip.t, abs(ip.data)*1e-6, 'k-','DisplayName','i_p (MA)');
  plot(ne.t, ne.data*1e-20, 'r-','DisplayName', 'n_e^{fir} (10*10^{19} m^{-3})')
  legend show;legend('location', 'best'); box on; 
  hold off;xticks([0 0.5 1 1.5 2.])

  ax2=subplot(2,2,2); hold on;
  plot(powers.nbi.t, powers.nbi.data*1e-6, 'k', 'DisplayName', 'NBI');
  try
    plot(powers.ec.t, powers.ec.data(:,10)*1e-6, 'r', 'DisplayName', 'EC');
  catch
    disp('No ECH in tree')
  end
  legend show; box on; 
  hold off;xticks([0 0.5 1 1.5 2.])
  ylabel('P (MW)')
  
  ax3=subplot(2,2,3); hold on;
  plot(z.t, z.data, 'k', 'DisplayName', 'z (m)');
  plot(delta.t, delta.data, 'r', 'DisplayName', '\delta');
  legend show; box on; 
  hold off;
  xlabel('t (s)'); xticks([0 0.5 1 1.5 2.])

  ax4=subplot(2,2,4); hold on;
  try
%       plot(b0.t, b0.data, 'k', 'DisplayName', 'Bphi')
    plot(theta_angle_launcher.t, theta_angle_launcher.data(1,:), 'k', 'DisplayName', 'L1 pol');
    plot(phi_angle_launcher.t, phi_angle_launcher.data(1,:), 'k--', 'DisplayName', 'L1 tor');
    plot(theta_angle_launcher.t, theta_angle_launcher.data(2,:), 'r', 'DisplayName', 'L2 pol');
    plot(phi_angle_launcher.t, phi_angle_launcher.data(2,:), 'r--', 'DisplayName', 'L2 tor');
  catch
      disp('No ECH angles in tree')
      title('No ECH in tree')
  end
  legend show; box on; 
  hold off;
  linkaxes([ax1,ax2,ax3, ax4],'x')
  xlim([0., 2.])
  xlabel('t (s)'); xticks([0 0.5 1 1.5 2.])
  
  %% save figure
  if save==1
      disp('saving file')
      saveas(gcf,sprintf('/home/vallar/missions/HLT11/W47/%i_summary.png', shotin(ishot)), 'png')
  end
end