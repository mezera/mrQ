function plotData(data, time, datafit, T1)

% Takes a 3D data set, where the third dimension is TIs
% and plots the different datapoints of a voxel selected with the mouse
    
% inspired by relaxPlotTimeSeries.m in mrvista
% (http://white.stanford.edu/software/)
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    
    
fig123 = figure(123); imshow(squeeze(data(:,:,1)), []);
movegui(fig123,'northwest')
button = 1;
while (button ==1)
   figure(123)
   [xx, yy, button] = ginput(1);
   yy = cast(floor(yy), 'int16');
   xx = cast(floor(xx), 'int16');
   fig = figure;
   movegui(fig,'northeast')
   x = 0:2500;
   if (nargin==2)
	   plot (time,squeeze(data(floor(yy), floor(xx), :)) , 'b+'); 
   end
   if (nargin>=3)
	   plot(time,squeeze(data(floor(yy), floor(xx), :)) , 'b+', linspace(min(time),max(time),20),squeeze(datafit(floor(yy), floor(xx),:)), 'r');
   end
   if(nargin==4)
       % Location is given from top left corner
       title (sprintf('Location X = %d, Y = %d, T1 = %g ms', cast(floor(yy), 'int16'), cast(floor(xx), 'int16'), T1(yy, xx)));
   end
   customFormat
   drawnow;
end