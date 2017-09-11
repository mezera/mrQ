function plot2Dhist(inX,inY,nBins,Xminmax,Yminmax,Xlable,Ylable,XTicks,YTicks)
% call to make a 2D histogram
%
%      Examle:
%
% inX=rand(1000,1);inY=rand(1000,1);nBins=30;Xminmax=[0 0.5];Yminmax=[0 0.5];Xlable='inX';Ylable='inY';XTicks=[0 1.25 0.4];YTicks=[0 1.25 0.4];
% plot2Dhist(inX,inY,nBins,Xminmax,Yminmax,Xlable,Ylable,XTicks,YTicks)


% Requires vistasoft, mrQ and knkutils
%
% AM/BW Vistaosft Team, 2013


if notDefined('nBins')
    nBins=100;
end

if notDefined('Xminmax')
    Xminmax=minmax(inX);
end

if notDefined('Yminmax')
    Yminmax=minmax(inY);
end




CV1=calccod(inX,inY,[],1);

[n,x,y] = mrQ_hist2d(inX,inY,nBins);maxN = ceil(max(n(:))/10)*10;


% figure; clf;

image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256).^3));
curve(1:2)=[min([Xminmax Yminmax]), max([Xminmax Yminmax])];
hold on;
h = plot(curve,curve,'--');
hold off;
axis([Xminmax Yminmax]);
axis square xy;
if ~notDefined('Xlable')
    xlabel(Xlable ,'FontSize',18);
end
if ~notDefined('Ylable')
    ylabel(Ylable ,'FontSize',18);
end

set(h,'LineWidth',2.5,'Color',[.4 .4 .4])%,'FontSize' ,'10');
set(gca,'FontSize',18)
if ~notDefined('XTicks')
    set(gca,'XTick',XTicks,'FontSize' ,14);
end
if ~notDefined('YTicks')
    set(gca,'YTick',YTicks,'FontSize' ,14),
end
% title(['R^2 = ' num2str(CV1)],'FontSize',18)


% colorbar ??