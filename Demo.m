load filename.txt;          % file with coordinates for the channel; retain file name
load loc_wells.txt;         % file with coordinates for any wells that may be present; retain file name; if there's no well leave content;
load property_SGSIM.txt;    % Point Bar Porosity Distribution Modeled Using SGSIM and ELLIPSIM

nx   =  150;                % Number of grid blocks along x-direction
ny   =  150;                % Number of grid blocks along y-direction
nz   =  25;                 % Number of grid blocks along z-direction
PBH  =  50;                 % Point Bar Height, ft
DTOP =  9973.753;           % depth to the top of the point bar reservoir, ft

[corner_point_grid,well_locs_locations,lac_x,lac_y,ihs_hor,ihs_vert] = point_bar_geomodel(filename,loc_wells,nx,ny,nz,PBH,DTOP,property_SGSIM);

% --- Dispay the GRID IN 2D
figure (3);
plot2Dgrid(lac_x,lac_y); 
title 'Curvilinear Grid (Aerial View)'

figure(4);
plot2Dgrid(ihs_hor,ihs_vert); 
title 'Curvilinear Grid (Sectional View)'



