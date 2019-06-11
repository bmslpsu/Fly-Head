function [] = GetBode()
%% GetBode:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
clear ; close all ; clc

T = 1;
Fs = 200;
t = (0:(1/Fs):T)';
tt = (0:(1/(Fs*100)):T)';
pp = (pi/4)*[0 1 2 3 4 5 6 7 8];
A = 1*ones(1,length(pp));
% ff = [0.5 1 2 3 4 5 6 7 8 9 10 11 12 13];
ff = 12*A;
nF = length(ff);
for kk = 1:nF
    U(:,kk) = A(kk)*cos(2*pi*ff(kk)*t + pp(kk));
%     U(:,kk) = spline(t,U(:,kk),tt);
    [Fv(:,kk), Mag(:,kk) , Phs(:,kk) , freq(:,kk)] = FFT(t,U(:,kk));
    rr(:,kk) = real(freq(:,kk));
    img(:,kk) = imag(freq(:,kk));
    [MAG(:,kk),PHASE(:,kk),fIdx(:,kk)] = Get_IO_Freq(Fv(:,kk),Mag(:,kk),Phs(:,kk),ff(kk));
    FREQ = freq(fIdx(:,kk),:);
 	R(:,kk) = real(FREQ(:,kk));
    I(:,kk) = imag(FREQ(:,kk));
end

pLabel = cellfun(@(x) num2str(rad2deg(x)),num2cell(pp),'UniformOutput',false);

figure (1) ; clf ; hold on
plot(t,U)
% plot(t,Y)

% figure (2) ; clf ; hold on
% h.circle = PlotCircle(0,0,1);
% h.circle.Color = 'k';
% h.circle.LineWidth = 1;
% h.cmplx = plot(rr,img,'*');
% % xlim([0 10])

Fig3 = figure (3) ; clf ; hold on
Fig3.Position = [100 100 600 600];
movegui(Fig3,'center')
ax = gca;
box on
grid on
axis equal
ax.XLim = 1.2*[-1 1];
ax.YLim = ax.XLim;
h.xorigin = plot((-1000:1000),0*(-1000:1000),'-g');
h.xorigin = plot(0*(-1000:1000),(-1000:1000),'-g');
h.circle = PlotCircle(0,0,1);
h.circle.Color = 'k';
h.circle.LineStyle = '-';
h.circle.LineWidth = 1;
h.circle = PlotCircle(0,0,0.5);
h.circle.Color = 'k';
h.circle.LineStyle = '-';
h.circle.LineWidth = 1;
for jj = 1:length(R)
    plotv([R(jj);I(jj)],'k')
    h.cmplx(jj) = plot(R(jj),I(jj),'.','MarkerSize',30);    
end
legend([h.cmplx],pLabel)
% plotv([rr,img]','-');





end
