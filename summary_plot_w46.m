function summary_plot_w46(shotin,tm,tM,doprint)
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
if ~exist('tm','var') || isempty(tm); tm=0;end;
if ~exist('tM','var') || isempty(tM); tM=2.;end;
if ~exist('doprint','var'); doprint=0;end;
if ~exist('freq_max','var'); freq_max=20;end;

if isempty(which('subtightplot'))
  tight_path = '/home/labit/matlab/mytcvlib';
  addpath(tight_path)
  cln = onCleanup(@() rmpath(tight_path));
end
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.1], [0.1 0.01], [0.06 0.06]);
if ~make_it_tight,  clear subplot;  end

% close all
% clear functions
set_defaults_matlab % to have color order from colos irrespective of matlab version
%set(0,'DefaultlineLineWidth',1); % do use width of 2 by default set in set_defaults_matlab
dc = get(gca,'colororder'); % Default colors
if abs(doprint)==1; close all; end % assume clean matlab better for printing mode
for ishot=1:length(shotin)
  shot = mdsopen(shotin(ishot));
  tL = mdsvalue('dim_of(\results::i_p)');
  % mdsclose;
  tM = nanmax(tM,nanmax(tL)+0.05);

  i=21;
  if abs(doprint)==1
    fig=figure('Visible','off','PaperPositionMode','auto');
  else
    fig=figure('Visible','on','PaperPositionMode','auto');
  end
  set(fig, 'PaperType', 'A4');
  font_weight='normal';
  if ~ldisplay
    % Optimize for fullHD resolution
    fontsize_norm = 20;
    pos = [1 1 1920 1080];
  else
    fontsize_norm = min(get(0,'DefaultAxesFontSize'),10);
    pos = [ 358 417 1085 690];
    pos = [ 358 417 1150 830];
  end
  fontsize_big  = fontsize_norm*3/2;
  set(fig,'position',pos);
  set(fig,'DefaultAxesFontsize',fontsize_norm,'DefaultTextFontSize',fontsize_norm);
  set(fig,'DefaultaxesFontWeight',font_weight,'DefaulttextFontWeight',font_weight);

  xte_auto_2000(shot); % Force to fill XTe nodes -- BL -- 2018.05.18

  [t,f] = get_gas_flux(shot,1,'D2'); f=f/1e20;
  cond = t>tm&t<tM;
  t(~cond)=[]; f(~cond)=[];
  [t2,f2] = get_gas_flux(shot,2,'N2'); f2=f2/1e20;
  cond = t2>tm&t2<tM;
  t2(~cond)=[]; f2(~cond)=[];
  [t3,f3] = get_gas_flux(shot,3,'N2'); f3=f3/1e20;
  cond = t3>tm&t3<tM;
  t3(~cond)=[]; f3(~cond)=[];

  mdsopen(shot);


  %% Powers
  subplot(4,2,5)
%  [tt,tmp] = tcvget('POHMSMOOTH');
%  cond = tt>tm & tt<tM;
%  plot(tt(cond),tmp(cond)/1e3,'Color',dc(1,:));
  hold all
  try
    main=tdi('\ATLAS::NBH.DATA.MAIN_ADC:DATA');
    cond = main.dim{1}>tm&main.dim{1}<tM;
  catch
    clear main
    main.dim{1}=linspace(0,2,101);
    main.data(:,37) = zeros(size(main.dim{1}));
    cond = main.dim{1}>tm&main.dim{1}<tM;
  end
  plot(main.dim{1}(cond),main.data(cond,37)*1000,'Color',dc(2,:)); %kW

  main=tdi('\results::toray:input:p_gyro[*,10]');
  cond = main.dim{1}>tm&main.dim{1}<tM;
  plot(main.dim{1}(cond),main.data(cond),'Color',dc(3,:)); %kW
  try
    main=tdi('\results::bolo:prad:power');
    cond = main.dim{1}>tm&main.dim{1}<tM;
  catch
    clear main
    main.dim{1}=linspace(0,2,101);
    main.data(:,37) = zeros(size(main.dim{1}));
    cond = main.dim{1}>tm&main.dim{1}<tM;
  end
  plot(main.dim{1}(cond),sum(main.data(cond,:),2),'Color',dc(4,:)); %kW
  axis tight;ylim = get(gca,'YLim'); set(gca,'YLim',[0,ylim(2)*1.1]);
  set(gca,'XTickLabel','','XGrid','On','YGrid','On');
  ylabel('Power (kW)')
  hl=legend('P_{ohm}','P_{NB}','P_{EC}','P_{rad}','Location','northwest'); set(hl,'Box','Off');


  %% Density and gas flux
  subplot(4,2,3)
%  main=tdi('\results::fir:n_average');
%  cond = main.dim{1}>tm&main.dim{1}<tM;
%  main.data = main.data/1e19;
%  ngf=gdat(shot,'ngf');
  a=plotyy(t,f,main.dim{1}(cond),main.data(cond)); %
  hold(a(1),'on')
  plot(a(1),t2,f2,'Color',dc(3,:),'Linewidth',2);plot(t3,f3,'Color',dc(4,:),'Linewidth',2);
  hold(a(2),'on')
%  plot(a(2),ngf.t,ngf.data*10);
  hl=legend(a(2),'n_{e,av19}','10*ngf','location','east'); set(hl,'Box','Off');
  %axis tight; ylim = get(a(1),'YLim'); ylim = ceil(ylim);
  %set(a(1),'YLim',[0 ylim(2)],'YTick',[0:2:ylim(2)],'YTickLabel',[0:2:ylim(2)]);
  %set(a(2),'YLim',[0 Inf])
  ylabel(a(1),'flux (x 1e20 mol/s) ') % ;ylabel(a(2),'n_{e,av} (x1e19 m-3)');
  set(a(1),'XTickLabel','','XGrid','On','YGrid','On');set(get(a(1),'Children'),'Linewidth',2);
  set(a(2),'XTickLabel','','XGrid','On','YGrid','On');set(get(a(2),'Children'),'Linewidth',2);
  aa=get(a(1),'Ytick');
  bb=get(a(2),'Ylim');
  set(a(2),'YLim',[0 bb(2)],'YTick',[0:bb(2)/(length(aa)-1):bb(2)]);

%   %% Stored energy
%   subplot(4,2,4)
% 
%   hold all
%   main=tdi('tcv_eq("total_energy","liuqe.m")');
%   try
%     cond = main.dim{1}>tm & main.dim{1}<tM;
%   catch
%     main.dim{1}=NaN; main.data=NaN; cond=true;
%   end
%   betan=gdat(shot,'betan');
%   aa=plotyy(main.dim{1}(cond),main.data(cond)/1e3,betan.t,betan.data);
%   set(aa(1),'Linewidth',2);
%   set(gca,'XTickLabel','','XGrid','On','Box','On','YGrid','On');
%   hl=legend('W_{MHD} (kJ)','\beta_N','Location','northwest'); set(hl,'Box','Off');


  %% Te(0)
  subplot(4,2,6)
  main=tdi('\results::thomson:profiles:auto:te_axis');
  try
    cond = isnan(main.data) | main.data<0;
    main.dim{1}(cond)=[]; main.data(cond)=[];
    cond = main.dim{1}>tm&main.dim{1}<tM;
  catch
    main.dim{1}=NaN; main.data=NaN; cond=true;
  end
  try
    xte=tdi('\results::te_x_a:foo');
    xte.data = xte.data(:,1);
    cx = xte.dim{1}>tm & xte.dim{1}<tM;
  catch
    clear xte
    xte.dim{1}=linspace(0,2,31);
    xte.data = NaN(size(xte.dim{1}));
    cx = xte.dim{1}>tm & xte.dim{1}<tM;
  end
  a=plot(xte.dim{1}(cx),xte.data(cx)/1e3,main.dim{1}(cond),main.data(cond)/1e3,'Linewidth',1); %keV
  set(a(2),'Linewidth',2);
  ylabel('T_{e,0} (keV)');
  hl=legend('XTE','TS','Location','southeast');set(hl,'Box','Off');
  try
    [zzmax,itzzmax] = nanmax(main.data(cond));
    it_xte=iround_os(xte.dim{1},main.dim{1}(cond(itzzmax)));
    set(gca,'XTickLabel','','XGrid','On','YGrid','On','YLim',[0.*nanmin(main.data(cond)/1e3) 1.1.*nanmax(xte.data(it_xte),zzmax)/1e3]);
  end
%   %% D-alpha
%   main=tdi('\pd:pd_001');
%   if sum(main.data>9.9)>1e3,
%     main =tdi('pd_calibrated(13)');
%   else
%     main =tdi('pd_calibrated( 1)');
%   end;
%   cond = main.dim{1}>tm&main.dim{1}<tM;
%   ax = subplot(4,2,7);
% %  q95=gdat(shot,'q95');
% %  a=plotyy(main.dim{1}(cond),main.data(cond),q95.t,q95.data); %
%   %plot(main.dim{1}(cond),main.data(cond),'Linewidth',1);
% % $$$ ylabel(ax,'D-alpha (cal.)')
% % $$$ xlabel(ax,'time (s)');
%   %set(gca,'YLim',[0 Inf],'XGrid','On','YGrid','On');
% %  set(a(1),'YLim',[0 Inf],'XGrid','On','YGrid','On');
% %  ylabel(a(1),'D-alpha (cal.)')
% %  ylabel(a(2),'q95')
% %  hhx=xlabel(a(1),'time (s)');
%   hold all
% % $$$ htemp=plot(q95.t,q95.data,'-');
%   hl=legend('D_\alpha','q_{95}','location','nw'); set(hl,'Box','Off');
% % $$$ text(tm,0,num2str(shot),'Fontsize',fontsize_big,'Fontweight','bold',...
% % $$$   'VerticalAlignment','Bottom','HorizontalAlignment','Left');
%   hhx_pos = get(hhx,'pos');
%   text(min(get(ax,'XLim')),hhx_pos(2),num2str(shot),'Fontsize',fontsize_big,'Fontweight','bold',...
%        'VerticalAlignment','Top','HorizontalAlignment','Left');

  %% BT and IP
  subplot(4,2,1)
  signB0 = mdsvalue('\pcs::mgams.data:if36fb');
  signIP = mdsvalue('\pcs::mgams.data:iohfb');
  bt = tdi('tcv_eq("BZERO")');
  cond = bt.dim{1}>tm&bt.dim{1}<tM;
  ip = tdi('tcv_eq("\magnetics::iplasma:trapeze")');
  a=plotyy(ip.dim{1},signIP*ip.data/1e3,bt.dim{1}(cond),signB0*bt.data(cond));
  %set(a(2),'YLim',[0 10])
  ylabel(a(1),'I_p (kA)');ylabel(a(2),'B_T (T)');
  set(a(1),'XTickLabel',''); set(get(a(1),'Children'),'Linewidth',2);
  set(a(2),'XTickLabel',''); set(get(a(2),'Children'),'Linewidth',1);
  set(a(1),'XGrid','on','YGrid','On')
  text(min(get(a(1),'XLim')),min(get(a(1),'YLim')),sprintf('I_p sign = %+1d, B_T sign = %+1d',signIP,signB0),...
       'Fontsize',fontsize_big,'Parent',a(1),...
       'VerticalAlignment','Bottom','HorizontalAlignment','Left');


  % % Soft-x and magnetics
  % try
  %     main = tdi('\atlas::dt100_northeast_001:channel_001');
  %     condm = main.dim{1}>tm&main.dim{1}<tM;
  % catch
  %     clear main
  %     main.dim{1}=linspace(0,2,31);
  %     main.data=zeros(size(main.dim{1}));
  %     condm= main.dim{1}>tm&main.dim{1}<tM;
  % end
  % try
  %     mag=tdi('\atlas::dt196_MHD_001:channel_083');
  %     [bb,aa] = butter(4,0.25);
  %     mag.data=filtfilt(bb,aa,mag.data);
  %     cond = mag.dim{1}>tm&mag.dim{1}<tM;
  % catch
  %     clear mag
  %     mag.dim{1}=linspace(0,2,31);
  %     mag.data=zeros(size(mag.dim{1}));
  %     cond= mag.dim{1}>tm&mag.dim{1}<tM;
  % end
  % subplot(4,2,8)
  % a=plotyy(main.dim{1}(condm),main.data(condm),mag.dim{1}(cond),mag.data(cond));
  % set(a(1),'YLim',[0 Inf],'XGrid','On','YGrid','On');set(a(2),'YLim',[-2 2])
  % ylabel(a(1),'SXR');ylabel(a(2),'magnetics');
  % xlabel('time (s)')


%   % Shaping
%   del = mdsvalue('tcv_eq("DELTA_EDGE","LIUQE.M")');
%   kap = mdsvalue('tcv_eq("KAPPA_EDGE","LIUQE.M")');
%   zmag = mdsvalue('tcv_eq("Z_AXIS","LIUQE.M")');
%   tL = mdsvalue('dim_of(tcv_eq("DELTA_EDGE","LIUQE.M"))');
%   cond = tL>tm & tL<tM;
%   subplot(4,2,2)
%   plot(tL(cond),zmag(cond),'Color',dc(1,:))
%   set(gca,'XTick',[],'XTickLabel','');
%   hold all
%   a=plotyy(tL(cond),del(cond),tL(cond),kap(cond));
% 
%   YTick2 = get(a(2),'YTick');
%   YTick = get(a(1),'YTick');
%   if numel(YTick)==2&&numel(YTick2)==2
%     YTick(3)=YTick(2); YTick(2)=mean([YTick(1),YTick(3)]);
%     YTick2(3)=YTick2(2); YTick2(2)=mean([YTick2(1),YTick2(3)]);
%   end
% 
%   set(a(2),'YTick',YTick2,'YTickLabel',num2str(YTick2'),'YColor',dc(3,:));
%   set(a(1),'YTick',YTick,'YTickLabel',num2str(YTick'),'YColor','k');
%   set(a(1),'XTickLabel','','XGrid','On','YGrid','On');
%   set(a(2),'XTickLabel','','XGrid','On','YGrid','On');set(get(a(2),'Children'),'Color',dc(3,:));
%   hl=legend('z_{mag}','d','k','Location','southeast');set(hl,'Box','Off');
% 
%   % Magnetics spectrogram
% 
%   if isempty(which('FMPloader')), MHDstartup;end
% 
%   if shot > 50921
%     array='POL-003';
%   else
%     array='POL';
%   end
% 
%   load /home/labit/matlab/mytcvlib/pol003.mat
%   mdsopen(shot); % Reopen shot
%   tmp = tdi('\results::z_axis');
%   % keyboard
%   zaxis = mean(tmp.data(tmp.dim{1}>=tm & tmp.dim{1}<=tM));
%   % mdsclose
% 
%   R0 = 0.88;
%   zpr = pol003.grid.x{2}; rpr = pol003.grid.x{1};
%   [~,ind] = min(abs(zpr-zaxis));
%   ds = abs(zpr(ind)-zaxis);
%   ind = find(abs(zpr-zaxis)<=2*ds&rpr>0.88,1,'first'); % Be sure that only one probe signal will be downloaded
%   %probes  = pol003.signalNames;
%   probes = str2num(pol003.signalNames{ind}(end-1:end));
% 
%   data3 = FMPloader('shot',shot,'array',array, 'calibrated',true,'tload',[tm tM],...
%           'probes',probes,'integrated',false,'blnRMoffset',true,'blnSat2nan',true,...
%           'blnBBtheta',false);
%   facq = 1/mean(diff(data3.tdi.dim{1}));
% 
%   [b,a] = butter(4,0.5);
%   for i=1:size(data3.tdi.data, 2)
%     data3.tdi.data(:,i) = filtfilt(b, a, data3.tdi.data(:,i));
%   end
%   data3.tdi.data = data3.tdi.data./repmat(std(data3.tdi.data,0,1),[size(data3.tdi.data,1),1]);
%   x = data3.tdi.data; t = data3.tdi.dim{1};
% 
% 
%   Nbuf = 2000;
%   clear XX
%   XX(:,:,1) = buffer(x(:,1), Nbuf, Nbuf-Nbuf/8, 'nodelay');
%   % XX(:,:,2) = buffer(x(:,end), Nbuf, Nbuf/2, 'nodelay');
% 
%   % XX = permute(XX,[1,3,2]);
%   tt = buffer(t, Nbuf, Nbuf-Nbuf/8, 'nodelay');
%   tt(:, end) = []; XX(:,  end) = []; tt = nanmean(tt, 1);
% 
%   clear Pxy
%   nfft=1024;
%   for i=1:size(XX,2)
%     [Pxy(:,i),f]=cpsd(XX(:,i),XX(:,i),hanning(nfft),0,nfft,facq/1e3);
%   end
% 
%   subplot(4,2,8)
%   cond=f<freq_max;
%   imagesc(tt,f(cond),log10(Pxy(cond,:))); set(gca,'YDir','normal');
%   set(gca,'ylim',[0 Inf]);
%   ylabel('f (kHz)');
%   xlabel('time (s)')
% 
% 
%   % a  =get(fig,'Children'); set(a(1),'Title',num2str(shot))
%   % mdsclose;
% 
% 
%   % Each subplot with the same x-axis
%   h=get(fig,'Children');
% 
%   if strcmp(version('-release'),'2014a')
%     set(h,'XLim',[tm tM])
%   else
%     for i=1:length(h)
%       if strcmp(h(i).Type,'axes')
%         set(h(i),'XLim',[tm tM]);
%       end
%     end
%   end
% 
%   % simple function to save snapshot figure to png
%   try
%     root_dir = '/Terra12/B/web/xtomoSVD';
%     if exist(root_dir,'dir'),
%       root_dir = fullfile(root_dir,[num2str(floor(shot./1e3)),'XXX']);
%     else
%       root_dir = fullfile('/tmp',getenv('USER'));
%     end
%     fn = fullfile(root_dir,['SUMMARY',num2str(shot),'.png']);
%     fn_local = fullfile('./',['SUMMARY',num2str(shot),'.png']);
% 
%     if ldisplay % Test, if there is a display
%       set(fig, 'PaperUnits', 'centimeters');
%       % set(fig, 'PaperPosition', [0 0 29 18]); %x_width=10cm y_width=15cm
%       set(gcf,'paperpositionmode','auto');
%       if abs(doprint)==1
%         saveas(fig,strrep(fn,'png','fig'),'fig');
%         saveas(fig,fn,'png');
%         %    fig2public('inp_folder',root_dir,'out_folder',root_dir,'other_formats',{'png'})
%         eval(['!rm -rf ',fullfile(root_dir,['SUMMARY',num2str(shot),'.fig'])]);
%         %     eval(['! rm -rf ',root_dir,'SUMMARY',num2str(shot),'.eps']);
%       elseif abs(doprint)==2
%         saveas(fig,fn_local,'png');
% % $$$         set(fig,'PaperOrientation','Portrait')
% % $$$         fn_local_eps = strrep(fn_local,'png','eps');
% % $$$         print(fig,'-loose','-depsc',fn_local_eps)
% % $$$         system(['LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH ','gs -dEPSCrop -q -dNOPAUSE -dBATCH', ...
% % $$$                 ' -r', num2str(sppi),' -sDEVICE=png16m',' -sOutputFile=', fn_local,' ', fn_local_eps]);
% % $$$         system(sprintf('convert %1$s -extent %2$dx%3$d %1$s',fn_local,1920,1080));
% % $$$         system(['rm ' fn_local_eps]);
%       end
%     else
%       % Borrowed from eq_plot
% 
%       % NODISPLAY mode: MATLAB  uses ghostscripts rather than
%       % MATLAB drivers to produce pngs. However within MATLAB
%       % Ghostscript ignores -r option.
% 
%       fn_eps = strrep(fn,'png','eps');
% 
%       set(fig,'PaperOrientation','Portrait')
%       if abs(doprint)==1
%         disp(['NODISPLAY: Print eps file and use gs ', ...
%               'to convert to png ...'])
%         print(fig,'-loose','-depsc',fn_eps)
%         % Convert to png with Full HD resolution
%         system(['LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH ',...
%                 'gs -dEPSCrop -q -dNOPAUSE -dBATCH', ...
%                 ' -r', num2str(sppi), ...
%                 ' -sDEVICE=png16m', ...
%                 ' -sOutputFile=', fn, ...
%                 ' ', fn_eps]);
% 
%         % Use ImageMagick to crop/pad with extra pixels for exact resolution
%         system(sprintf('convert %1$s -extent %2$dx%3$d %1$s',fn,1920,1080));
% 
%         disp([' ... and remove temporary file ', fn_eps]);
%         system(['rm ' fn_eps]);
%       elseif abs(doprint)==2
%         fn_local_eps = strrep(fn_local,'png','eps');
%         print(fig,'-loose','-depsc',fn_local_eps)
%         system(['LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH ','gs -dEPSCrop -q -dNOPAUSE -dBATCH', ...
%                 ' -r', num2str(sppi),' -sDEVICE=png16m',' -sOutputFile=', fn_local,' ', fn_local_eps]);
%         system(sprintf('convert %1$s -extent %2$dx%3$d %1$s',fn_local,1920,1080));
%         system(['rm ' fn_local_eps]);
%       end
%     end
% 
%     if shot == 0 && abs(doprint)==1,
%       unix(sprintf('ln -fv /Terra12/B/web/xtomoSVD/%02dXXX/SUMMARY%05d.png /Terra12/B/web/xtomoSVD/SUMMARYLS.png',floor(shot/1000),shot));
%     end
%     fprintf('\n        =======================\n        === shot %d done ===\n        =======================\n\n',shot);
% 
%     if doprint < 0; close(fig); end
% 
%   catch ME
%     if abs(doprint)==1; disp('.. Saving figure at /Terra12/B/web/xtomoSVD/ failed'); end
%     disp(ME.getReport);
%   end


end

if ~strcmp(getenv('USER'),'labit')
  % (on cleanup) % rmpath /home/labit/matlab/mytcvlib
end
