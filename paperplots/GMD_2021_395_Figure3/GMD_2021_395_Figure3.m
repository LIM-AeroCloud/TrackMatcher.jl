clear all; close all;

%% Import data from text file, intercepts from Tesche et al. (2016)
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["lon", "lat", "alt", "timediff", "timediffhh", "day0night1"];
opts.VariableTypes = ["double", "double", "double", "double", "datetime", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "timediffhh", "InputFormat", "HH:mm:ss");

% Import the data
tbl = readtable("GMD_2021_395_Figure3_interceptstesche.txt", opts);

%% Convert to output type
lon1 = tbl.lon;
lat1 = tbl.lat;
alt1 = tbl.alt;
timediff = tbl.timediff;
timediffhh = tbl.timediffhh;
day0night1 = tbl.day0night1;

%% Clear temporary variables
clear opts tbl


%% Import data from text file, intercepts from TrackMatcher
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["lon", "lat", "alt"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("GMD_2021_395_Figure3_interceptstrackmatcher.txt", opts);

%% Convert to output type
lon2 = tbl.lon;
lat2 = tbl.lat;
alt2 = tbl.alt;

%% Clear temporary variables
clear opts tbl

lon1 = lon1 +180;
lon2 = lon2 +180;

lon_binID =  1 + floor( lon1 ./ 1.0 ) ;
lat_binID = 180 + floor( lat1 ./ 0.5 ) ;
binID     = [lon_binID, lat_binID] ;
nBin      = 360 ; 
onesVector = ones( size(lon1,1), 1 ) ;
counts1  = accumarray( binID, onesVector, [nBin, nBin] ) ;

lon_binID =  1 + floor( lon2 ./ 1.0 ) ;
lat_binID = 180 + floor( lat2 ./ 0.5 ) ;
binID     = [lon_binID, lat_binID] ;
nBin      = 360 ; 
onesVector = ones( size(lon2,1), 1 ) ;
counts2  = accumarray( binID, onesVector, [nBin, nBin] ) ;

difference = counts2 - counts1;

% plot the data
load coastlines;
newmap = parula;                    
newmap(1,:) = [1 1 1];
bins = [0:1:360];
c=[1 10/3 10 100/3 100 1000/3 1000];
limit = max(max(counts1));
  
s = get(0, 'ScreenSize'); 
f = figure('Position', [0 0 s(3) s(4)]);
subplot(2,1,1) 
  worldmap([16 48],[-160 -110])
  h2 = surfm(bins.*0.5-90,bins.*1.0-180,counts2'); 
  set(h2,'edgecolor','none'); 
  hold on;
  plotm(coastlat,coastlon,'-k');
  hold on;
  caxis([1 150])
%   colorbar;
  cb = colorbar();
  cb.Ruler.Scale = 'log';
  cb.Ruler.MinorTick = 'on';
  colormap(gca,newmap);
  title('Occurrence rate of intercepts, 2.5 h delay, 0.5 by 1.0 grid cell, TrackMatcher')

subplot(2,1,2) 
  worldmap([16 48],[-160 -110])
  h2 = surfm(bins.*0.5-90,bins.*1.0-180,difference'); 
  set(h2,'edgecolor','none'); % remove edges
  hold on;
  plotm(coastlat,coastlon,'-k');
  hold on;
  caxis([-170 170])
  cb = colorbar();
  cb.Ruler.Scale = 'lin';
  cb.Ruler.MinorTick = 'on';
  colormap(gca,redblue);
  title('Occurrence rate of intercepts, 2.5 h delay, 0.5 by 1.0 grid cell, difference TrackMatcher - Tesche et al. (2016)')

set(f,'PaperOrientation','landscape');
% print(f, '-dpdf', '-painters', '-fillpage', 'GMD_2021_395_Figure3.pdf');
% print(f, '-depsc2', '-tiff', '-r300', '-painters', 'GMD_2021_395_Figure3.eps'); 
% set(f,'PaperOrientation','portrait');
% print(f, '-dpng', '-r150', 'GMD_2021_395_Figure3.png');