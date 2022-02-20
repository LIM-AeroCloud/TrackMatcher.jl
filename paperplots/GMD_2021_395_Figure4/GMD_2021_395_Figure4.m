clear all; close all;

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["id", "lat", "lon", "alt", "tdiff", "tflight", "tsat", "atmos_state"];
opts.VariableTypes = ["double", "double", "double", "double", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["tdiff", "tflight", "tsat"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["tdiff", "tflight", "tsat", "atmos_state"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "id", "TrimNonNumeric", true);
opts = setvaropts(opts, "id", "ThousandsSeparator", ",");

% Import the data
X2012021011 = readtable("GMD-2021-395-Figure4.csv", opts);

%% Convert to output type
lon = X2012021011.lon;
lat = X2012021011.lat;
alt = X2012021011.alt;
timediff = X2012021011.tdiff;
cloud = X2012021011.atmos_state;

f = find(cloud == 'ci');

%% Clear temporary variables
clear opts tbl

lon = lon +180;

% transform to grid, 2.5 degree by 1.25 degree in this case
lon_binID =  1 + floor( lon ./ 2.5 ) ;
lat_binID = 72 + floor( lat ./ 1.25 ) ;
binID     = [lon_binID, lat_binID] ;
nBin      = 144 ; 
onesVector = ones( size(lon,1), 1 ) ;
counts  = accumarray( binID, onesVector, [nBin, nBin] ) ;

lon_binID2 =  1 + floor( lon(f) ./ 2.5 ) ;
lat_binID2 = 72 + floor( lat(f) ./ 1.25 ) ;
binID2     = [lon_binID2, lat_binID2] ;
nBin      = 144 ; 
onesVector2 = ones( size(lon(f),1), 1 ) ;
counts2  = accumarray( binID2, onesVector2, [nBin, nBin] ) ;



% plot the data
load coastlines;
% replace 0 with white in color map
newmap = parula;                    
newmap(1,:) = [1 1 1];
bins = [0:1:144];

s = get(0, 'ScreenSize'); 
f = figure('Position', [0 0 s(3) s(4)]);
subplot(2,1,1) 
  worldmap([-90 90],[-180 180])
  h2 = surfm(bins.*1.25-90,bins.*2.5-180,counts'); 
  set(h2,'edgecolor','none'); 
  hold on;
  plotm(coastlat,coastlon,'-k');
  hold on;
  caxis([1 200])
  cb = colorbar();
  cb.Ruler.Scale = 'log';
  cb.Ruler.MinorTick = 'on';
  colormap(newmap);
  title('Occurrence rate of intercepts above 5 km height, 0.5 h delay, 0.5 by 1.0 grid cell')

subplot(2,1,2)   
  worldmap([-90 90],[-180 180])
  h2 = surfm(bins.*1.25-90,bins.*2.5-180,counts2'); 
  set(h2,'edgecolor','none'); 
  hold on;
  plotm(coastlat,coastlon,'-k');
  hold on;
  caxis([1 30])
  cb = colorbar();
  cb.Ruler.Scale = 'log';
  cb.Ruler.MinorTick = 'on';
  colormap(newmap);
  title('Occurrence rate of intercepts above 5 km height, 0.5 h delay, 0.5 by 1.0 grid cell, within cirrus')

set(f,'PaperOrientation','landscape');
print(f, '-dpdf', '-painters', '-fillpage', 'GMD-2021-395-Figure4');
print(f, '-depsc2', '-tiff', '-r300', '-painters', 'GMD-2021-395-Figure4'); 
% set(f,'PaperOrientation','portrait');
% print(f, '-dpng', '-r150', 'GMD-2021-395-Figure4.png');