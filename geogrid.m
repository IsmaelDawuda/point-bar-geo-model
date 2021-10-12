function [corner_point_grid,well_locs_locations,lac_x,lac_y,ihs_hor,ihs_vert] = geogrid(filename,loc_wells,nx,ny,nz,PBH,DTOP)

% ----- CODE BEGINS -----

loc_wells = loc_wells';

ChannelPath = filename; % 
channelx    =  ChannelPath(:,1);
channely    =  ChannelPath(:,2);

x_y   =  sortrows([channelx';channely'].',1).';  % % channel path sorted in the direction of channel progression
xout  =  x_y(1,:);
yout  =  x_y(2,:);



channel_x = xout;  
channel_y = yout; 


ch_points      =  [channel_x;channel_y]; 


% --- ROTATION OF POINTS ----%

rotation_angle = input("Enter a rotation angle to align points with the global X-Y plane, Enter 0 if points are already aligned: "); %angle is  +ve for anticlockwise rotation and -ve for clockwise

% ROTATION MATRIX
R = [cosd(rotation_angle) -sind(rotation_angle); sind(rotation_angle) cosd(rotation_angle)];
% Rotate your point(s)
trans_pts = zeros(2,numel(ch_points(1,:)));

for i = 1:numel(ch_points(1,:))
trans_pts(:,i) =  R*ch_points(:,i);
end


channel_path_x  =  trans_pts(1,:);
channel_path_y  =  trans_pts(2,:); 


[pks, pklocs]    = findpeaks(channel_path_y);                               % Original Peaks
[troughs,trlocs] = findpeaks(-channel_path_y);                           % Original Troughs

pktr    =  [pklocs(:), pks(:); trlocs(:) -troughs(:)];          % Location, Value Matrix
pktrs   =  sortrows(pktr);                                     % Location, Value Matrix Sorted By Location
rr      =  pktrs(:,1);

pktrsx  =  channel_path_x(rr);
pktrsy  =  pktrs(:,2);
pktrsy  =  pktrsy';


cr_tr_y =  [channel_path_y(1),pktrsy,channel_path_y(end)];

cr_tr_x =  [channel_path_x(1),pktrsx,channel_path_x(end)];

 
  slp   =  zeros(1,numel(cr_tr_x)-1);
for i   =  1:numel(cr_tr_x)-1
slp(i)  =  (cr_tr_y(i)-cr_tr_y(i+1))./(cr_tr_x(i)-cr_tr_x(i+1));
end


ii = numel(slp)-1;
if numel(slp) == 2
pbar_num = 1;
else
pbar_num = input(['Which of the ' num2str(ii) ' concave parts should be used? the number should be <= number of concave parts: ']);
end

conc_num =(numel(cr_tr_x)-2);
assert(( pbar_num <= conc_num ),'Out of Bounds---Number Should not be more than the number of concave parts of the channel')


colm = find(channel_path_x >= cr_tr_x(pbar_num) & channel_path_x <= cr_tr_x(pbar_num+2));

x = channel_path_x(colm);
y = channel_path_y(colm); 


if slp(pbar_num)<slp(pbar_num+1)
CONCAVE_DIR = 2;
elseif slp(pbar_num)>slp(pbar_num+1)
CONCAVE_DIR = 1;
end



N       = 20000;
points  = [x;y];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               


%----sine curve conditioned to a cubic spline --------% 
A     =  fnplt(cscvn(points));

newx  =  A(1,:); newy  =  A(2,:); 

trans_well_locs = zeros(2,numel(loc_wells(1,:)));

for i = 1:numel(loc_wells(1,:))
   trans_well_locs(:,i) =  R*loc_wells(:,i);
end
loc_wells = trans_well_locs;
loc_wells_x = trans_well_locs(1,:);
loc_wells_y = trans_well_locs(2,:);

col_m = find(loc_wells_x > x(1) & loc_wells_x < x(end));


w_x = [x x(end) x(1)];  % closing points to form a polygon
w_y = [y y(end) y(1)]; % closing points to form a polygon
%plot(w_x,w_y)
xqw = loc_wells_x(col_m);   %%x-cordinates of well locations in and out of the polygon
yqw = loc_wells_y(col_m);   %% y-cordinates of well locations in and out of the polygon
%hold on; scatter(xqw,yqw);
in  = inpolygon(xqw,yqw,w_x,w_y);  


well_locs_x = xqw(in); % x-cordinates of well locations within the concave part of the meander
well_locs_y = yqw(in); % y-cordinates of well locations within the concave part of the meander
well_locs   = [well_locs_x ; well_locs_y];

% plot(x,y); hold on;
% scatter(well_locs_x,well_locs_y,'r','filled') % points inside

 % ----- INJECTION OR PRODUCING WELL LOCATIONS -----------% 
 
num_of_inj_wells  =     numel(well_locs(1,:)); 
well_xcord        =     well_locs_x;
well_ycord        =     well_locs_y;

w   =  input('What is the angular displacement of this channel? :  ');

% ---CHECK WHETHER ANGLE SPECIFIED IS VALID FOR SINE GENERATION FUNCTION----

minVal     =  0;
maxVal     =  89.5;
assert(( minVal<w) && ( w<maxVal ),'Out of Bounds---Angle  Specified Should Satisfy the Condition: 0 < w < 89.5')

 
 myd = zeros(1,numel(newx));  myl = zeros(1,numel(newx));

for i = 1:numel(newx)
    
if CONCAVE_DIR == 2 && newy(2)< newy(1)
    
w       =   -w;
xstart  =  newx(1);
ystart  =  newy(1);
yend    =  ystart;
ymid    =  min(y);
f       =  find(newy==ymid);
xmid    =  newx(f);  xmid = min(xmid);
dfr     =  (xmid - xstart);      
xend    =  xstart + 2*dfr;

elseif CONCAVE_DIR ==2 && newy(end)> newy(end-1)
w             =   -w;
ymid          =  min(newy);
g             =  find(newy==ymid);
xmid          =  newx(g); xmid  = max(xmid);
yend          =   newy(end);
xend          =  newx(end);
ystart        =  yend;
xstart        =  xend - 2*abs((xend - xmid));



elseif CONCAVE_DIR ==1 && newy(end)< newy(end-1) 
 
ymid          =  max(newy);
ff            =  find(newy==ymid);
xmid          =  newx(ff); xmid  = max(xmid);
yend          =  min(newy);
rf            =  find(newy==yend);
xend          =  newx(end);
ystart        =  yend;

xstart        =  xend - 2*(xend - xmid);



elseif CONCAVE_DIR ==1 && newy(2)> newy(1)
 
xstart       =  newx(1);
ystart       =  newy(1);
yend         =  ystart;
ymid         =  max(y);
f            =  find(newy==max(newy));
xmid         =  newx(f);  xmid = min(xmid);
dfr          =  abs(xmid - xstart);      
xend         =  xstart + 2*dfr;

end

    thet    =  (w/180)*pi;

    
    if i==1
        myd  = 0;
        myl  = 0; 
       
    else
    
          myd(i)    =  sqrt(((newx(i)-newx(i-1)).^2) + ((newy(i)-newy(i-1)).^2));
          myl(i)    =  sum(myd(1:i));  
   
          
    end
end


    
myL         =  max(myl);
    
xcord_fun2  =  @(myl) cos(thet*sin(2*pi*myl/myL));
ycord_fun2  =  @(myl) sin(thet*sin(2*pi*myl/myL));

my_xcord2 =zeros( 1,numel(newx)); my_ycord2 = my_xcord2;
for i = 1:numel(newx)
  
          my_xcord2(i) = xstart + integral(xcord_fun2,0,myl(i));          % ---convert to cartesian coordinates ----
   
          my_ycord2(i) = ystart + integral(ycord_fun2,0,myl(i));          % ---convert to cartesian coordinates ----
          
          
end


lambda  = sqrt((my_xcord2(1) - my_xcord2(end)).^2  + (my_ycord2(1) - my_ycord2(end)).^2) ;  % wavelength

   
% ------------ FIND FOCAL POINT ----------


% ------------FIND SLOPES BTN EACH MEANDER AND THE FOCUS, CALC THE CORDINATES AT THE TOE OF THE CHANNELS-------

 
 %%lamda  =   10.9*(W.^1.01);     %(Leopold and Wolman, 1960) units is in meters
 %lam     =   (lambda/10.9);
 lam     =  (lambda/10.9);
 W       =   nthroot(lam,1.01);
 
Wo           =   W/1.5;                 %(Leopold and Wolman, 1960), this unit is in meters
 
conv_factor  =  3.28084;    
 d       =   (conv_factor)*Wo;     %  unit converted to feets
 % REFERENCE:  Leopold, L. B., and M. G. Wolman, 1960, River meanders: Geological Society of America Bulletin, v. 71, p. 769 – 794

 s2 = zeros(1,numel(my_ycord2)); mytheta2 = s2; x2_toe = s2; y2_toe=s2;
 
for i = 1:numel(my_ycord2)
if CONCAVE_DIR ==2 
   y_foc   =   max(my_ycord2);
else
   y_foc   =   min(my_ycord2);
end

   x_foc        =   median(my_xcord2);   


   s2(i)        =  (my_ycord2(i) - y_foc)./ (my_xcord2(i) - x_foc);   % slope for  current position of channel  
   mytheta2(i)  =   atand(s2(i));
   
   
  if s2(i)    >=0 && CONCAVE_DIR ==2
  
  x2_toe(i)  =  my_xcord2(i) - d.*cosd(mytheta2(i));  
  y2_toe(i)  =  my_ycord2(i) - d.*sind(mytheta2(i));
  
  
  elseif s2(i)  < 0 && CONCAVE_DIR ==2
  
  
  x2_toe(i)  =  my_xcord2(i) + d.*cosd(mytheta2(i));  
  y2_toe(i)  =  my_ycord2(i) + d.*sind(mytheta2(i));
  
  
  
  elseif s2(i)    >=0 && CONCAVE_DIR ==1
  
  x2_toe(i)  =  my_xcord2(i) + d.*cosd(mytheta2(i));  
  y2_toe(i)  =  my_ycord2(i) + d.*sind(mytheta2(i));
  
  
  elseif s2(i)  < 0 && CONCAVE_DIR ==1
  
  
  x2_toe(i)  =  my_xcord2(i) - d.*cosd(mytheta2(i));  
  y2_toe(i)  =  my_ycord2(i) - d.*sind(mytheta2(i));
  
  
  end
end




    


 % intial channel position 
nc  =  100;
xc1  =  zeros(nc,numel(my_xcord2)); 
yc1  =  zeros(nc,numel(my_ycord2));


for i = 1:numel(my_xcord2)
   xc1(:,i) = linspace(my_xcord2(i),x_foc,nc);
   
   yc1(:,i) = linspace(my_ycord2(i),y_foc,nc);
end


my_xcord1  =  xc1(10,:);  
my_ycord1  =  yc1(80,:);
 
 

 
 
 
 x1_toe = zeros(1,numel(my_ycord2)); y1_toe = x1_toe;
for i = 1:numel(my_ycord2)
      
  if s2(i)    >=0 &&CONCAVE_DIR ==2
  
  x1_toe(i)  =  my_xcord1(i) - d.*cosd(mytheta2(i));  
  y1_toe(i)  =  my_ycord1(i) - d*sind(mytheta2(i));
  
  
  elseif s2(i) <0 &&CONCAVE_DIR ==2
  
  
  x1_toe(i)  =  my_xcord1(i) + d.*cosd(mytheta2(i));  
  y1_toe(i)  =  my_ycord1(i) + d.*sind(mytheta2(i));
  
  
  
  
  elseif s2(i)    >=0 &&CONCAVE_DIR ==1
  
  x1_toe(i)  =  my_xcord1(i) + d.*cosd(mytheta2(i));  
  y1_toe(i)  =  my_ycord1(i) + d*sind(mytheta2(i));
  
  
  elseif s2(i) <0 && CONCAVE_DIR ==1
  
  
  x1_toe(i)  =  my_xcord1(i) - d.*cosd(mytheta2(i));  
  y1_toe(i)  =  my_ycord1(i) - d.*sind(mytheta2(i));
  
  
  end
end


% ---remove duplicates  ---

my_xcord1  =  my_xcord1(3:3:end);  
my_ycord1  =  my_ycord1(3:3:end);
my_xcord2  =  my_xcord2(3:3:end);  
my_ycord2  =  my_ycord2(3:3:end);
x2_toe     =  x2_toe(3:3:end);
y2_toe     =  y2_toe(3:3:end);

x1_toe     =  x1_toe(3:3:end);
y1_toe     =  y1_toe(3:3:end);
    


% -----------ACCOUNT FOR THE RELATIVE LOCATION OF THE WELLS TO THE SHIFTED CHANNEL   -----------%

% calculate distances
D_HOR     =  abs(xstart - xend);   % horizontal distance between the start and end

D_HOR_FNL =  abs(my_xcord2(1)- my_xcord2(end)); % horizontal distance between the start and end

D_VERT    = abs(.5*(ystart+yend) - ymid);   

D_VERT_FNL    = abs(.5*(my_ycord2(1)+my_ycord2(end)) - min(my_ycord2));  


D_VERT_FNLL    = abs(max(my_ycord2)- min(my_ycord2));  


%find the closest value of xcoordinate to the well location;
d_hor = zeros(1,num_of_inj_wells);  d_vert = d_hor; ratio_hor = d_hor; ratio_vert = d_hor;
well_ycord_fnl =d_hor; well_xcord_fnl = d_hor;
for i =1:num_of_inj_wells
d_hor(i)            =   abs(well_xcord(i) - xstart);

d_vert(i)           =   abs(well_ycord(i) - .5*(ystart +yend));


ratio_hor(i)        =   d_hor(i)/D_HOR;

if CONCAVE_DIR      ==  2
ratio_vert(i)       =   d_vert(i)/D_VERT;



elseif CONCAVE_DIR ==  1
    
D_VERT_FNL          =  D_VERT_FNLL;
ratio_vert(i)       =   (1- (d_vert(i)/D_VERT));

end

well_ycord_fnl(i)   =  max(my_ycord2)  - ratio_vert(i)*(D_VERT_FNL);
well_xcord_fnl(i)   =  min(my_xcord2)    + (ratio_hor(i)*(D_HOR_FNL));


end

x1 = my_xcord2; y1 =  my_ycord2;

x2 = my_xcord1; y2 =  my_ycord1;

 

xx1   =  linspace(x1(1),x1(end),N);   % query points along x for one side of the channel
yy1   =  spline(x1,y1,xx1);           % corresponding y value for every query  point along x for one side of the channel

xx1_toe   =  linspace(x2_toe(1),x2_toe(end),N);  
yy1_toe   =  spline(x2_toe,y2_toe,xx1_toe);           



xx2   =  linspace(x2(1),x2(end),N);  % query points along x for the other side of the channel
yy2   =  spline(x2,y2,xx2);          % corresponding y value for every query  point along x for the other side of the channel

xx2_toe    =  linspace(x1_toe(1),x1_toe(end),N);  
yy2_toe    =  spline(x1_toe,y1_toe,xx2_toe);          




% ---- calc the true length along the channel path   for one side of the channel
l1 = zeros(1,numel(x1)); d1 = zeros(1,numel(x1));  D1 = zeros(1,numel(x1));  d1_toe = l1; l1_toe=l1;
for i = 1:numel(xx1)
   
    if i==1
         d1(i)     =  0;
         l1(i)     =  0;
         
         d1_toe(i) =  0;
         l1_toe(i) =  0;

    
    
    else
         
         d1(i)       =  sqrt(((xx1(i)-xx1(i-1)).^2) + ((yy1(i)-yy1(i-1)).^2));
         
         d1_toe(i)   =  sqrt(((xx1_toe(i)-xx1_toe(i-1)).^2) + ((yy1_toe(i)-yy1_toe(i-1)).^2));
         
         
         l1(i)       =  sum(d1(1:i)); % the total distance covered at any point along the channel path 
         
         l1_toe(i)   =  sum(d1_toe(1:i)); 
         
        Channel_length_1 = max(l1); % this is the total distance covered at the end of the channel path

        Channel_length_1_toe = max(l1_toe); 

    end
end



l2 = zeros(1,numel(x2)); d2 = zeros(1,numel(x2));    D2 = zeros(1,numel(x2)); d2_toe =l2; l2_toe=l2;
% ---- calc the true length along the channel path   for the side of the channel
for i = 1:numel(xx2)
   
    if i==1
         d2(i)=0;
         l2(i)=0;
         
         d2_toe(i)=0;
         l2_toe(i)=0;
         
    else
        
         d2(i) = sqrt(((xx2(i)-xx2(i-1)).^2) + ((yy2(i)-yy2(i-1)).^2));
         d2_toe(i) = sqrt(((xx2_toe(i)-xx2_toe(i-1)).^2) + ((yy2_toe(i)-yy2_toe(i-1)).^2));
         
         l2(i) = sum(d2(1:i)); % the total distance covered at any point along the channel path 
         l2_toe(i) = sum(d2_toe(1:i)); % the total distance covered at any point along the channel path 
         
         
          Channel_length_2 = max(l2); % this is the total distance covered at the end of the channel path
          Channel_length_2_toe = max(l2_toe);
         
    
    
    end
end


blks_alongcurve           =    nx;                         %   number of grid blocks along curve curve
N_B                       =    ny;                           % blocks btn the two gridded surfaces
blks_along_sigmoid        =    nz;   
points_along_curve        =    blks_alongcurve + 1;           % the no of grid nodes along each curve
blksize_curve1            =    Channel_length_1/blks_alongcurve;
blksize_curve1_toe        =    Channel_length_1_toe/blks_alongcurve;
blksize_curve2            =    Channel_length_2/blks_alongcurve;
blksize_curve2_toe        =    Channel_length_2_toe/blks_alongcurve;

cumdist1 = l1; cumdist2 = l2;

cumdist1_toe = l1_toe; cumdist2_toe = l2_toe;

D1_toe = zeros(points_along_curve,1); D2_toe = D1_toe;
for i  =  1:points_along_curve
D1(i)  =   (i-1)*Channel_length_1/blksize_curve1 ; % total distance traveled at every node along meander for the lower part
D1_toe(i)  =   (i-1)*Channel_length_1_toe/blksize_curve1_toe ; 

D2(i)      =   (i-1)*Channel_length_2/blksize_curve2 ; % total distance traveled at every node along meander for the upper part
D2_toe(i)  =   (i-1)*Channel_length_2_toe/blksize_curve2_toe ;

end

discretise1       =   cumdist1(1):blksize_curve1:cumdist1(end);                % divisions into blocks
discretise1_toe   =   cumdist1_toe(1):blksize_curve1_toe:cumdist1_toe(end);    % divisions into blocks
discretise2       =   cumdist2(1):blksize_curve2:cumdist2(end);                % divisions into blocks
discretise2_toe   =   cumdist2_toe(1):blksize_curve2_toe:cumdist2_toe(end);    % divisions into blocks

% ----compare discretise with total distance covered and choose the nearest number in the total distance () -----
 TMP1        =   bsxfun(@(x1,y1) abs(x1-y1), discretise1(:), reshape(cumdist1,1,[]));
 TMP1_toe    =   bsxfun(@(x1_toe,y1_toe) abs(x1_toe-y1_toe), discretise1_toe(:), reshape(cumdist1_toe,1,[]));
 TMP2        =   bsxfun(@(x2,y2) abs(x2-y2), discretise2(:), reshape(cumdist2,1,[]));
 TMP2_toe    =   bsxfun(@(x2_toe,y2_toe) abs(x2_toe-y2_toe), discretise2_toe(:), reshape(cumdist2_toe,1,[]));
 
 [~, idxcumdist1]           =   min(TMP1,[],2);
 [D1_toe, idxcumdist1_toe]   =   min(TMP1_toe,[],2);
[D2, idxcumdist2]            =   min(TMP2,[],2);
[D2_toe, idxcumdist2_toe]    =   min(TMP2_toe,[],2);
Result1                      =   cumdist1(idxcumdist1);
Result1_toe                  =   cumdist1_toe(idxcumdist1_toe);
Result2                      =   cumdist2(idxcumdist2);
Result2_toe                  =   cumdist2_toe(idxcumdist2_toe);


 
%---find where the distance along curve matches the total distance for descretisation 
search1 = zeros(1,points_along_curve);search2 =   zeros(1,points_along_curve);
search1_toe = search1; search2_toe = search1;
for i = 1:numel(Result1)
search1(i)  = find(cumdist1 ==  Result1(i));

search1_toe(i)  = find(cumdist1_toe ==  Result1_toe(i));

search2(i)  = find(cumdist2 ==  Result2(i));

search2_toe(i)  = find(cumdist2_toe ==  Result2_toe(i));
end
%  ---extract the x y coordinates for every division 
%for each distance traveled we have the x y coordinates

xcord1    =  xx1(search1); % x cordinate for curve 1
xcord2    =  xx2(search2); % x cordinate for curve 2
ycord1    =  yy1(search1); % y cordinate for curve 1
ycord2    =  yy2(search2); % y cordinate for curve 2


xcord1_toe    =  xx1_toe(search1_toe); 
xcord2_toe    =  xx2_toe(search2_toe); 
ycord1_toe    =  yy1_toe(search1_toe);
ycord2_toe    =  yy2_toe(search2_toe);



% --- Specify the number of blocks btn the two surfaces ----
             % the number of block btn the two surfaces
np          =   blks_along_sigmoid +  1;  % no of grid nodes btn the blks 

%-----Get the distances and angles   between each pair of the extracted coordinates ----
dist = zeros(1,numel(xcord1));  m = zeros(1,numel(xcord1)); ang = zeros(1,numel(xcord1)); blks_size_pepn = zeros(1,numel(xcord1));
mydist = dist; mym = dist; myang = dist; blks_size_btn_surf =dist;
dist_btm = dist;  Mm = dist;  Aang = dist;  blks_size_pepn_btm =dist;
for i    = 1:numel(xcord1)
mydist(i)  = sqrt((xcord1_toe(i) - xcord2(i)).^2 + (ycord1_toe(i) - ycord2(i)).^2) ; 
mym(i)     = (ycord1_toe(i) - ycord2(i)) / (xcord1_toe(i) - xcord2(i)) ;
myang(i)   = atand(mym(i));
blks_size_btn_surf(i)  =  mydist(i)/blks_along_sigmoid;    % block size btn the two surfaces

% top surf ----
dist(i)  = sqrt((xcord1_toe(i) - xcord1(i)).^2 + (ycord1_toe(i) - ycord1(i)).^2) ; 
m(i)     = (ycord1_toe(i) - ycord1(i)) / (xcord1_toe(i) - xcord1(i)) ;
ang(i)   = atand(m(i));


blks_size_pepn(i)  =  dist(i)/blks_along_sigmoid;    % block size btn the two surfaces



% ---  bottom surf ----


dist_btm(i)  = sqrt((xcord2_toe(i) - xcord2(i)).^2 + (ycord2_toe(i) - ycord2(i)).^2) ; 
Mm(i)     = (ycord2_toe(i) - ycord2(i)) / (xcord2_toe(i) - xcord2(i)) ;
Aang(i)   = atand(Mm(i));


blks_size_pepn_btm(i)  =  dist_btm(i)/blks_along_sigmoid;    % block size btn the two surfaces

end


 myx_ihs_top_surf = zeros(numel(xcord1),np); myy_ihs_top_surf = zeros(numel(xcord1),np);
 myx = myx_ihs_top_surf; myy = myx_ihs_top_surf; myx_btm = myx_ihs_top_surf; myy_btm = myx_ihs_top_surf; 
 for i   = 1:numel(xcord1)
for j   = 1:np
    
if m(i) > 0 && CONCAVE_DIR==2

 myx_ihs_top_surf(i,j)  = xcord1_toe(i) + (j-1)*blks_size_pepn(i)*cosd(ang(i));
 myy_ihs_top_surf(i,j)  = ycord1_toe(i) + (j-1)*blks_size_pepn(i)*sind(ang(i));


elseif m(i) < 0 && CONCAVE_DIR==2

myx_ihs_top_surf(i,j)  = xcord1_toe(i) - (j-1)*blks_size_pepn(i)*cosd(ang(i));
myy_ihs_top_surf(i,j)  = ycord1_toe(i) - (j-1)*blks_size_pepn(i)*sind(ang(i));

end


if m(i) > 0 && CONCAVE_DIR==1

 myx_ihs_top_surf(i,j)  = xcord1_toe(i) - (j-1)*blks_size_pepn(i)*cosd(ang(i));
 myy_ihs_top_surf(i,j)  = ycord1_toe(i) - (j-1)*blks_size_pepn(i)*sind(ang(i));


elseif m(i) < 0 && CONCAVE_DIR==1

myx_ihs_top_surf(i,j)  = xcord1_toe(i) + (j-1)*blks_size_pepn(i)*cosd(ang(i));
myy_ihs_top_surf(i,j)  = ycord1_toe(i) + (j-1)*blks_size_pepn(i)*sind(ang(i));

end



if myang(i)>0   && CONCAVE_DIR==2
 myx(i,j)  = xcord1_toe(i) - (j-1)*blks_size_btn_surf(i)*cosd(myang(i));
 myy(i,j)  = ycord1_toe(i) - (j-1)*blks_size_btn_surf(i)*sind(myang(i));
 
elseif myang(i)< 0 && CONCAVE_DIR==2
myx(i,j)  = xcord1_toe(i) +(j-1)*blks_size_btn_surf(i)*cosd(myang(i));
myy(i,j)  = ycord1_toe(i) + (j-1)*blks_size_btn_surf(i)*sind(myang(i));

end


if myang(i)>0   && CONCAVE_DIR==1
 myx(i,j)  = xcord1_toe(i) + (j-1)*blks_size_btn_surf(i)*cosd(myang(i));
 myy(i,j)  = ycord1_toe(i) + (j-1)*blks_size_btn_surf(i)*sind(myang(i));
 
elseif myang(i)< 0 && CONCAVE_DIR==1
myx(i,j)  = xcord1_toe(i) -(j-1)*blks_size_btn_surf(i)*cosd(myang(i));
myy(i,j)  = ycord1_toe(i) - (j-1)*blks_size_btn_surf(i)*sind(myang(i));

end

if Aang(i) >  0 && CONCAVE_DIR==2   
    
myx_btm(i,j)  = xcord2_toe(i) + (j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
 myy_btm(i,j)  = ycord2_toe(i) + (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
elseif Aang(i) <  0 && CONCAVE_DIR==2   
myx_btm(i,j)  = xcord2_toe(i) -(j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
myy_btm(i,j)  = ycord2_toe(i) - (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
    
    

end

if Aang(i) > 0 && CONCAVE_DIR==1  
    
myx_btm(i,j)  = xcord2_toe(i) - (j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
 myy_btm(i,j)  = ycord2_toe(i) - (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
elseif Aang(i) <0 && CONCAVE_DIR ==1
myx_btm(i,j)  = xcord2_toe(i) +(j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
myy_btm(i,j)  = ycord2_toe(i) + (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
    
    

end



end
end






for i   = 1:numel(xcord1)
for j   = 1:np
if m(i) > 0 && CONCAVE_DIR==2 

 myx_ihs_top_surf(i,j)  = xcord1_toe(i) + (j-1)*blks_size_pepn(i)*cosd(ang(i));
 myy_ihs_top_surf(i,j)  = ycord1_toe(i) + (j-1)*blks_size_pepn(i)*sind(ang(i));


elseif m(i) < 0 && CONCAVE_DIR==2

myx_ihs_top_surf(i,j)  = xcord1_toe(i) - (j-1)*blks_size_pepn(i)*cosd(ang(i));
myy_ihs_top_surf(i,j)  = ycord1_toe(i) - (j-1)*blks_size_pepn(i)*sind(ang(i));

end


if m(i) > 0 && CONCAVE_DIR==1 

 myx_ihs_top_surf(i,j)  = xcord1_toe(i) - (j-1)*blks_size_pepn(i)*cosd(ang(i));
 myy_ihs_top_surf(i,j)  = ycord1_toe(i) - (j-1)*blks_size_pepn(i)*sind(ang(i));


elseif m(i) < 0 && CONCAVE_DIR==1

myx_ihs_top_surf(i,j)  = xcord1_toe(i) + (j-1)*blks_size_pepn(i)*cosd(ang(i));
myy_ihs_top_surf(i,j)  = ycord1_toe(i) + (j-1)*blks_size_pepn(i)*sind(ang(i));

end




if myang(i)>0 && CONCAVE_DIR ==2
 myx(i,j)  = xcord1_toe(i) + (j-1)*blks_size_btn_surf(i)*cosd(myang(i));
 myy(i,j)  = ycord1_toe(i) + (j-1)*blks_size_btn_surf(i)*sind(myang(i));
elseif myang(i)<0 && CONCAVE_DIR ==2
myx(i,j)  = xcord1_toe(i) -(j-1)*blks_size_btn_surf(i)*cosd(myang(i));
myy(i,j)  = ycord1_toe(i) - (j-1)*blks_size_btn_surf(i)*sind(myang(i));

end

if myang(i)>0 && CONCAVE_DIR ==1
 myx(i,j)  = xcord1_toe(i) - (j-1)*blks_size_btn_surf(i)*cosd(myang(i));
 myy(i,j)  = ycord1_toe(i) - (j-1)*blks_size_btn_surf(i)*sind(myang(i));
elseif myang(i)<0 && CONCAVE_DIR ==1
myx(i,j)  = xcord1_toe(i) +(j-1)*blks_size_btn_surf(i)*cosd(myang(i));
myy(i,j)  = ycord1_toe(i) + (j-1)*blks_size_btn_surf(i)*sind(myang(i));

end



if Aang(i) > 0   && CONCAVE_DIR ==2
    
myx_btm(i,j)  = xcord2_toe(i) + (j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
 myy_btm(i,j)  = ycord2_toe(i) + (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
elseif Aang(i) < 0 && CONCAVE_DIR ==2
myx_btm(i,j)  = xcord2_toe(i) -(j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
myy_btm(i,j)  = ycord2_toe(i) - (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
    
    

end


if Aang(i) > 0   && CONCAVE_DIR ==1
    
myx_btm(i,j)  = xcord2_toe(i) - (j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
 myy_btm(i,j)  = ycord2_toe(i) - (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
elseif Aang(i) < 0 && CONCAVE_DIR ==1
myx_btm(i,j)  = xcord2_toe(i) +(j-1)*blks_size_pepn_btm(i)*cosd(Aang(i));
myy_btm(i,j)  = ycord2_toe(i) + (j-1)*blks_size_pepn_btm(i)*sind(Aang(i));
    
                                                             

end


end
end




% ----GET THE  SIGMOIDS----
%PBH            =   50;                   % POINT BAR HEGHT
P              =   50;
xsigm1         =   linspace(2,6,P);   xsigm2       =   linspace(3,10,P); 

Sigm_incx_top  =  linspace(x2_toe(1),x1(1),P);
Sigm_decx_top  =  linspace(x2_toe(end),x1(end),P);

Sigm_incx_top_bt  =  linspace(x1_toe(1),x2(1),P);
Sigm_decx_top_bt  =  linspace(x1_toe(end),x2(end),P);
a_fnl        =   2;                         
a            =   3;
c1           =   4;  c2 = 6;
Z_inc        =   PBH*(sigmf(xsigm1,[a c1]));
Z_dec        =   PBH*sigmf(xsigm1,[-a c1]);

Z_inc_bt     =   PBH*sigmf(xsigm2,[a_fnl c2]); 
Z_dec_bt     =   PBH*sigmf(xsigm2,[-a_fnl c2]); 

% generate grid nodes along the various nodes

sigm_incx_top   =  linspace(x2_toe(1),x1(1),N);  % query points along x for the other side of the channel
z_inc           =  spline(Sigm_incx_top,Z_inc,sigm_incx_top);    


sigm_decx_top   =  linspace(x2_toe(end),x1(end),N);  % query points along x for the other side of the channel
z_dec           =  spline(Sigm_decx_top,Z_dec,sigm_decx_top);          % corresponding y value for every query  point along x for the other side of the channel

  l3 = zeros(1,numel(x1)); d3 = zeros(1,numel(x1));  D1 = zeros(1,numel(x1)); d4 = zeros(1,N); l4 = d4;
for i = 1:N
   
    if i==1
         d3(i)=0;
         l3(i)=0;
         d4(i)=0;
         l4(i)=0;
    else

         d3(i) = sqrt(((sigm_incx_top(i)-sigm_incx_top(i-1)).^2) + ((z_inc(i)-z_inc(i-1)).^2));
         
         l3(i) = sum(d3(1:i)); 
                    
         d4(i) = sqrt(((sigm_decx_top(i)-sigm_decx_top(i-1)).^2) + ((z_dec(i)-z_dec(i-1)).^2));
         
         l4(i) = sum(d4(1:i)); % the total distance covered at any point along the channel path           
     end
end

sigm_length_2 = max(l4);
sigm_length_1 = max(l3);
points_alongsigmoid         =    blks_along_sigmoid + 1;           % the no of grid nodes along each curve
blksize_sigm_inc            =    sigm_length_1/blks_along_sigmoid;
blksize_sigm_dec            =    sigm_length_2/blks_along_sigmoid;

cumdist3 = l3; cumdist4 = l4;

D3 = zeros(1,points_alongsigmoid);  D4 = zeros(1,points_alongsigmoid);
for i  =  1:points_alongsigmoid
D3(i)  =   (i-1)*blksize_sigm_inc  ; % total distance traveled at every node along meander for the lower part
D4(i)  =   (i-1)*blksize_sigm_dec ; % total distance traveled at every node along meander for the upper part
end

discretise3  = cumdist3(1):blksize_sigm_inc:cumdist3(end); % divisions into blocks
discretise4  = cumdist4(1):blksize_sigm_dec:cumdist4(end); % divisions into blocks

% ----compare discretise with total distance covered and choose the nearest number in the total distance () -----
 TMP3 = bsxfun(@(x3,y3) abs(x3-y3), discretise3(:), reshape(cumdist3,1,[]));
 TMP4 = bsxfun(@(x4,y4) abs(x4-y4), discretise4(:), reshape(cumdist4,1,[]));
[D3, idxcumdist3] = min(TMP3,[],2);
[D4, idxcumdist4] = min(TMP4,[],2);
Result3 = cumdist3(idxcumdist3);
Result4 = cumdist4(idxcumdist4);

% TFDiffLessThen3 = D1 < 3;
search3 = zeros(1,points_alongsigmoid);search4 =   zeros(1,points_alongsigmoid); 
for i = 1:points_alongsigmoid
search3(i) = find(cumdist3==Result3(i));
search4(i)  = find(cumdist4==Result4(i));
end
%  ---extract the x y coordinates for every division 
%for each distance traveled we have the x y coordinates

myz_inc    =  z_inc(search3); % x cordinate for curve 1
xsigminc   =  sigm_incx_top(search3);
myz_dec    =  z_dec(search4);
xsigmdec   =  sigm_decx_top(search4);
 

myz = zeros(points_along_curve,np);

% -- Assign these to the x y gridded nodes;

for i = 2:points_along_curve
  for j = 1:np

   if myx(i,j)>myx(i-1,j) % it's increasing
      myz(i,j) = myz_inc(j); 
  else
      myz(i,j) = myz_dec(j);
  end
  end
end
myz(1,:) = myz_inc;



%SIGMOIDS%

% % generate grid nodes along the various nodes

sigm_incx_top_bt   =  linspace(x1_toe(1),x2(1),N);  % query points along x for the other side of the channel
z_inc_bt           =  spline(Sigm_incx_top_bt,Z_inc_bt,sigm_incx_top_bt);    

sigm_decx_top_bt   =  linspace(x1_toe(end),x2(end),N);  % query points along x for the other side of the channel
z_dec_bt           =  spline(Sigm_decx_top_bt,Z_dec_bt,sigm_decx_top_bt);          % corresponding y value for every query  point along x for the other side of the channel


% ---- calc the true length along the channel path   for one side of the channel

d4_bt = zeros(1,N); l4_bt = d4_bt; l3_bt = zeros(1,N); d3_bt = zeros(1,N);
for i = 1:N
   
    if i==1
         d3_bt(i)=0;
         l3_bt(i)=0;
         d4_bt(i)=0;
         l4_bt(i)=0;
    else

         d3_bt(i) = sqrt(((sigm_incx_top_bt(i)-sigm_incx_top_bt(i-1)).^2) + ((z_inc_bt(i)-z_inc_bt(i-1)).^2));
         
         l3_bt(i) = sum(d3_bt(1:i)); 
               
         d4_bt(i) = sqrt(((sigm_decx_top_bt(i)-sigm_decx_top_bt(i-1)).^2) + ((z_dec_bt(i)-z_dec_bt(i-1)).^2));
         
         l4_bt(i) = sum(d4_bt(1:i)); % the total distance covered at any point along the channel path     
     end
end

sigm_length_2_bt = max(l4_bt);
sigm_length_1_bt = max(l3_bt);


       %   number of grid blocks along curve curves
points_alongsigmoid         =    blks_along_sigmoid + 1;           % the no of grid nodes along each curve
blksize_sigm_inc_bt         =    sigm_length_1_bt/blks_along_sigmoid;
blksize_sigm_dec_bt         =    sigm_length_2_bt/blks_along_sigmoid;

cumdist3_bt = l3_bt; cumdist4_bt = l4_bt;

D3_bt = zeros(points_alongsigmoid);  D4_bt = D3_bt; 
for i  =  1:points_alongsigmoid
D3_bt(i)  =   (i-1)*blksize_sigm_inc_bt  ; % total distance traveled at every node along meander for the lower part
D4_bt(i)  =   (i-1)*blksize_sigm_dec_bt ; % total distance traveled at every node along meander for the upper part
end


discretise3_bt  = cumdist3_bt(1):blksize_sigm_inc_bt:cumdist3_bt(end); % divisions into blocks
discretise4_bt  = cumdist4_bt(1):blksize_sigm_dec_bt:cumdist4_bt(end); % divisions into blocks

% ----compare discretise with total distance covered and choose the nearest number in the total distance () -----
TMP3_bt = bsxfun(@(x3_bt,y3_bt) abs(x3_bt-y3_bt), discretise3_bt(:), reshape(cumdist3_bt,1,[]));
TMP4_bt = bsxfun(@(x4_bt,y4_bt) abs(x4_bt-y4_bt), discretise4_bt(:), reshape(cumdist4_bt,1,[]));
[D3_bt, idxcumdist3_bt] = min(TMP3_bt,[],2);
[D4_bt, idxcumdist4_bt] = min(TMP4_bt,[],2);
Result3_bt = cumdist3_bt(idxcumdist3_bt);
Result4_bt = cumdist4_bt(idxcumdist4_bt);

% TFDiffLessThen3 = D1 < 3;
search3_bt = zeros(1,points_alongsigmoid);search4_bt =   zeros(1,points_alongsigmoid); 
for i = 1:points_alongsigmoid
search3_bt(i) = find(cumdist3_bt==Result3_bt(i));
search4_bt(i)  = find(cumdist4_bt==Result4_bt(i));
end
%  ---extract the x y coordinates for every division 
%for each distance traveled we have the x y coordinates

myz_inc_bt    =  z_inc_bt(search3_bt); % x cordinate for curve 1
xsigminc_bt   =  sigm_incx_top_bt(search3_bt);
myz_dec_bt    =  z_dec_bt(search4_bt);
xsigmdec_bt   =  sigm_decx_top_bt(search4_bt);
 
% -- Assign these to the x y gridded nodes;
myz_bt = zeros(points_along_curve,np);
for i = 2:points_along_curve
  for j = 1:np

   if myx_btm(i,j)>myx_btm(i-1,j) % it's increasing
      myz_bt(i,j) = myz_inc_bt(j); 
  else
      myz_bt(i,j) = myz_dec_bt(j);
  end
  end
end
 myz_bt(1,:) = myz_inc_bt;
 % ----Get the x y z coordinates in between
% calculate the coordinates btn the toe and top;

N_P   =  N_B + 1;                        % points btn the two gridded surfaces

DPN   =  zeros(points_along_curve,N_P,points_alongsigmoid);
DIST  =  zeros(points_along_curve,N_P,points_alongsigmoid);

for i = 1:points_along_curve
    for j = 1:N_P
        for k = 1:points_alongsigmoid
        
DIST(i,j,k)  = sqrt((myx_ihs_top_surf(i,k) - myx_btm(i,k)).^2 + (myy_ihs_top_surf(i,k) - myy_btm(i,k)).^2) ; 

DPN(i,j,k)  = DIST(i,j,k)/N_B;    % block size btn the two surfaces
        end
    end
end

% calculate the coordinates btn the toe and top;
finalx = zeros(points_along_curve,N_P,points_alongsigmoid);    
finaly = zeros(points_along_curve,N_P,points_alongsigmoid);
finalz = zeros(points_along_curve,N_P,points_alongsigmoid);

M = zeros(points_along_curve,points_alongsigmoid ); my_ang=M;
for i   = 1:points_along_curve
for j   = 1:N_P
for k   = 1:points_alongsigmoid  
  
M(i,k)        =   (myy_ihs_top_surf(i,k)- myy_btm(i,k))/(myx_ihs_top_surf(i,k)- myx_btm(i,k));
my_ang(i,k)   =   atand(M(i,k));    
if M(i,k) > 0 &&  CONCAVE_DIR==2
finalx(i,j,k)  =  myx_ihs_top_surf(i,k)  +    (j-1)*DPN(i,j,k)*cosd(my_ang(i,k));
finaly(i,j,k)  =  myy_ihs_top_surf(i,k)  +    (j-1)*DPN(i,j,k)*sind(my_ang(i,k));
elseif M(i,k) < 0 && CONCAVE_DIR==2
finalx(i,j,k)  =  myx_ihs_top_surf(i,k)  -    (j-1)*DPN(i,j,k)*cosd(my_ang(i,k));
finaly(i,j,k)  =  myy_ihs_top_surf(i,k)  -    (j-1)*DPN(i,j,k)*sind(my_ang(i,k));

end

if M(i,k) > 0 && CONCAVE_DIR==1
finalx(i,j,k)  =  myx_ihs_top_surf(i,k)  -    (j-1)*DPN(i,j,k)*cosd(my_ang(i,k));
finaly(i,j,k)  =  myy_ihs_top_surf(i,k)  -    (j-1)*DPN(i,j,k)*sind(my_ang(i,k));
elseif M(i,k) < 0 && CONCAVE_DIR==1
finalx(i,j,k)  =  myx_ihs_top_surf(i,k)  +    (j-1)*DPN(i,j,k)*cosd(my_ang(i,k));
finaly(i,j,k)  =  myy_ihs_top_surf(i,k)  +    (j-1)*DPN(i,j,k)*sind(my_ang(i,k));

end

finalz(i,j,k)  = myz(i,k) - (j-1)*(abs(myz(i,k) - myz_bt(i,k))/N_B);

end
end
end



lac_x = finalx(:,:,end);   lac_x = reshape(lac_x,nx+1,ny+1);
lac_y = finaly(:,:,end);   lac_y = reshape(lac_y,nx+1,ny+1);

ihs_hor = finaly(end,:,:);   ihs_hor = reshape(ihs_hor,ny+1,nz+1);
ihs_vert = finalz(end,:,:);   ihs_vert = reshape(ihs_vert,ny+1,nz+1);


XXX = reshape(finalx,[],1); XXX =  (XXX);
YYY = reshape(finaly,[],1); YYY =  (YYY);
ZZZ = reshape(finalz,[],1); ZZZ =  (ZZZ);

AAA = [XXX YYY ((flipud(ZZZ)) )];

wz  =PBH*(ones(1,num_of_inj_wells)); % just to show the location of wells at the top

myfinalx = reshape(finalx,N_P*points_along_curve,points_alongsigmoid); myfinaly = reshape(finaly,N_P*points_along_curve,points_alongsigmoid);
myfinalz = reshape(finalz,N_P*points_along_curve,points_alongsigmoid);


% %-------IMPLEMENT GRID TRANSFORMATION SCHEME%-------------
% % alph    =      0; %general orientation of the channel with respect to the horizontal plane
% % bet     =      90 - alph;  %general orientation of the channel perpendicular to the downstream direction
% %gam     =        %general orientation of the sigmoid
% % determine distance along acc surfs accounting for the orientations
if Channel_length_1>Channel_length_2 
DD      =     Channel_length_1;
else
    DD   = Channel_length_2;
end

% % determine distance along digmoid accounting for the orientations
dsy =  finaly(floor((N_P)/2),floor((N_P)/2),:);
dsz =  finalz(floor((N_P)/2),floor((N_P)/2),:);
dsy = dsy(:); 
dsz = dsz(:);

SS = zeros(1,points_alongsigmoid);
for i = 1:points_alongsigmoid
    if i ==1
        SS(i) =0;
    
    else
SS(i) = sqrt(((dsy(i)-dsy(i-1)).^2) + ((dsz(i)-dsz(i-1)).^2));
    end
end
ss = sum(SS(1:end));
% determine distance between acc surfs accounting for the orientations
dbn     =      max(dist);

[X, Y, Z] = ndgrid((linspace(0,DD,points_along_curve)), (linspace(0,dbn,N_P)), (linspace(0,ss,points_alongsigmoid)));
X = X+xcord1(1);
Y = Y+ycord1(1);

XX  =  reshape(X,[],1);  YY  = reshape(Y,[],1);  ZZ  =  reshape(Z,[],1);  % cordinates are converted  into column vectors

AA = [XX YY ((flipud(ZZ)) + 3040)];


% ----  END ------



% ---PLOTTING THE RECTILINEAR GRID ----%
xx_beg = X(:,1,:);
yy_beg = Y(:,1,:);
zz_beg = Z(:,1,:);

xx_end = X(:,end,:);
yy_end = Y(:,end,:);
zz_end = Z(:,end,:);

XX_beg = X(1,:,:);
YY_beg = Y(1,:,:);
ZZ_beg = Z(1,:,:);

XX_end = X(end,:,:);
YY_end = Y(end,:,:);
ZZ_end = Z(end,:,:);

xx_up  = X(:,:,end);
yy_up  = Y(:,:,end);
zz_up  = Z(:,:,end);

xx_down= X(:,:,1);
yy_down= Y(:,:,1);
zz_down =Z(:,:,1);


aspx = 1.5; aspy=1;aspz = .5;


% ----  END ------




%CONVERT CORNER POINT GRID TO CENTER GRID FOR BOTH CURVILOINEAR AND RECTILINEAR GRID 

xcoord_center_curv = zeros(points_along_curve-1,N_P -1,points_alongsigmoid -1);
ycoord_center_curv = zeros(points_along_curve-1,N_P -1,points_alongsigmoid -1);
zcoord_center_curv = zeros(points_along_curve-1,N_P -1,points_alongsigmoid -1);

xcoord_center_rect = zeros(points_along_curve-1,N_P -1,points_alongsigmoid -1);
ycoord_center_rect = zeros(points_along_curve-1,N_P -1,points_alongsigmoid -1);
zcoord_center_rect = zeros(points_along_curve-1,N_P -1,points_alongsigmoid -1);
for i   = 1:points_along_curve-1
for j   = 1:N_P -1
for k   = 1:points_alongsigmoid -1  
xcoord_center_curv(i,j,k) = (1/8)*(  finalx(i,j,k) +  finalx(i+1,j,k) + finalx(i+1,j+1,k)  + finalx(i,j+1,k) +  finalx(i+1,j,k+1)  +  finalx(i,j,k+1) +  finalx(i,j+1,k+1) +  finalx(i+1,j+1,k+1) );

ycoord_center_curv(i,j,k) = (1/8)*(  finaly(i,j,k) +  finaly(i+1,j,k) + finaly(i+1,j+1,k)  + finaly(i,j+1,k) +  finaly(i+1,j,k+1)  +  finaly(i,j,k+1)  +  finaly(i,j+1,k+1) +  finaly(i+1,j+1,k+1) );

zcoord_center_curv(i,j,k) = (1/8)*(  finalz(i,j,k) +  finalz(i+1,j,k) + finalz(i+1,j+1,k)  + finalz(i,j+1,k) +  finalz(i+1,j,k+1)  +  finalz(i,j,k+1)  +  finalz(i,j+1,k+1) +  finalz(i+1,j+1,k+1) );


xcoord_center_rect(i,j,k) = (1/8)*(  X(i,j,k) +  X(i+1,j,k) + X(i+1,j+1,k)  + X(i,j+1,k) +  X(i+1,j,k+1)  +  X(i,j,k+1)  +  X(i,j+1,k+1) +  X(i+1,j+1,k+1) );

ycoord_center_rect(i,j,k) = (1/8)*(  Y(i,j,k) +  Y(i+1,j,k) + Y(i+1,j+1,k)  + Y(i,j+1,k) +  Y(i+1,j,k+1)  +  Y(i,j,k+1) +  Y(i,j+1,k+1) +  Y(i+1,j+1,k+1) );

zcoord_center_rect(i,j,k) = (1/8)*(  Z(i,j,k) +  Z(i+1,j,k) + Z(i+1,j+1,k)  + Z(i,j+1,k) +  Z(i+1,j,k+1)  +  Z(i,j,k+1) +  Z(i,j+1,k+1) +  Z(i+1,j+1,k+1) );

end
end
end


xcmg                =   reshape(finalx,[],1);
ycmg                =   reshape(finaly,[],1);
zcmg                =   reshape(finalz,[],1);
coordinates_cmg     =   [xcmg ycmg zcmg]; 
my_coordinates     =   [xcmg ycmg  (flipud(zcmg))];  

%  ---- BACK TRANSFORM COORDINATES
rect_corner_grid = [XX YY ZZ];
rect_center_grid = [xcoord_center_rect(:) ycoord_center_rect(:) flipud(zcoord_center_rect(:)) ];
curv_center_grid = [xcoord_center_curv(:) ycoord_center_curv(:) flipud(zcoord_center_curv(:)) ];



curv_center_grid = curv_center_grid'; 
rect_center_grid = rect_center_grid';
back_trans_rect_cen  = zeros(3,numel(rect_center_grid(1,:)));  bk_trans_curv_cen = back_trans_rect_cen;
for i = 1:numel(curv_center_grid(1,:))
  rotation_angle = -rotation_angle;
  back_trans_rect_cen(1:2,i) =  R*(rect_center_grid(1:2,i));
  bk_trans_curv_cen(1:2,i) = R*(curv_center_grid(1:2,i));
end
back_trans_rect_cen(3,:) = rect_center_grid(3,:);   % THE Z COORDINATES HAVE NOT CHANGED
bk_trans_curv_cen(3,:) = curv_center_grid(3,:);   % THE Z COORDINATES HAVE NOT CHANGED

my_rect_center = back_trans_rect_cen';
my_curv_center = bk_trans_curv_cen';

%my_curv_center(:,3)  = flipud(my_curv_center(:,3));
save('rect_center_grid.txt','my_rect_center','-ascii')

save('curv_center_grid.txt','my_curv_center','-ascii')



my_coordinates = my_coordinates'; 
rect_corner_grid = rect_corner_grid';
back_trans_pts  = zeros(3,numel(my_coordinates(1,:)));  bk_trans_rect_corner = back_trans_pts;
for i = 1:numel(my_coordinates(1,:))
  rotation_angle = -rotation_angle;
  back_trans_pts(1:2,i) =  R*(my_coordinates(1:2,i));
  bk_trans_rect_corner(1:2,i) = R*(rect_corner_grid(1:2,i));
end

back_trans_pts(3,:) = my_coordinates(3,:);   % THE Z COORDINATES HAVE NOT CHANGED
bk_trans_rect_corner(3,:) = rect_corner_grid(3,:);

corner_point_grid = back_trans_pts';       % final grid generated in corner point format;
%corner_point_grid(:,3) = flipud(corner_point_grid(:,3)); 
%corner_point_grid = (1/conv_factor)*corner_point_grid;  % final converted back to feets;
my_rect_corner = bk_trans_rect_corner';
save('curv_corner_grid.txt','corner_point_grid','-ascii');




save('rect_corner_grid.txt','my_rect_corner','-ascii');
%save('rect_center_grid.txt','rect_center_grid','-ascii');





% BACK TRANSFORM THE WELL LOCATIONS TOO
my_new_well_locs = [well_xcord_fnl;well_ycord_fnl;wz];
back_trans_well_locs  = zeros(3,numel(my_new_well_locs(1,:)));
for i = 1:numel(my_new_well_locs(1,:))
  rotation_angle = -rotation_angle;
  back_trans_well_locs(1:2,i) =  R*(my_new_well_locs(1:2,i));
end


well_locs_locations = back_trans_well_locs;  % well locations converted back to feets

save('well_locs_locations.txt','back_trans_well_locs','-ascii');


Resx    =  reshape(finalx, points_along_curve*N_P,points_alongsigmoid);
Resy    =  reshape(finaly, points_along_curve*N_P,points_alongsigmoid);
Resz    =  reshape(finalz, points_along_curve*N_P,points_alongsigmoid);
 


 
figure(1)
scatter(filename(:,1),filename(:,2),'b','filled')  % initial points supplied
aa = 1:numel(filename(:,1));aa=aa'; b = num2str(aa); c = cellstr(b);
 dx = 20; dy =20; % displacement so the text does not overlay the data points, it can be adjusted depending on the data given
%text(channel_x+dx4, channel_y + dy4, c);
text(filename(:,1)+dx, filename(:,2)+dy, c);

hold on
plot(filename(:,1),filename(:,2),'b')  % initial points supplied
legend ('Channel Nodes','Unsorted Meander Path'); title 'Raw Data'
hold off

figure(2)
scatter(channel_x,channel_y,'b','filled');
a = 1:numel(channel_x);a=a'; b = num2str(a); c = cellstr(b);
 dxx = 20; dyy =20; % displacement so the text does not overlay the data points, it can be adjusted depending on the data given
%text(channel_x+dx4, channel_y + dy4, c);
text(channel_x+dxx, channel_y+dyy, c);

legend ('Sorted Channel Nodes'); title 'channel nodes sorted in the direction of channel progression'



BB     =  fnplt(cscvn(trans_pts));
sorted_ch_path_x  =  BB(1,:); sorted_ch_path_y  =  BB(2,:); 
figure(3)
for  jj =1:ii
colm = find(sorted_ch_path_x >= cr_tr_x(jj) & sorted_ch_path_x <= cr_tr_x(jj+2));

Xpol = sorted_ch_path_x(colm);
x_pol = [Xpol Xpol(end) Xpol(1)];

Ypol = sorted_ch_path_y(colm); 
y_pol = [Ypol Ypol(end) Ypol(1)];

npbar = 50;    % number of possible point bar well locations to simulate
x_q = x_pol(1) + (max(x_pol)-min(x_pol)) .* rand(npbar,1);
y_q = min(y_pol) + (max(y_pol)-min(y_pol)).*rand(npbar,1);
in  = inpolygon(x_q,y_q,x_pol,y_pol);

scatter(channel_x,channel_y,'b','filled');
a = 1:numel(channel_x);a=a'; b = num2str(a); c = cellstr(b);
 dx = 20; dy =20; % displacement so the text does not overlay the data points, it can be adjusted depending on the data given
%text(channel_x+dx4, channel_y + dy4, c);
t = text(channel_x+dx, channel_y+dx, c);

% hold on
% for kk =1:numel(a)
% t = text(channel_x+dx, channel_y+dx, c); t(kk).FontSize =14;
% end

hold on
plot(Xpol,Ypol,'b') 
hold on
scatter(x_q(in),y_q(in),'r','filled') % points inside
title 'meander path going through channel nodes and bending to accomodate possible point Bar well locations on the concave side of the meander'
legend ('Channel Nodes','Sorted Channel Path','Possible Point Bar Well Locations') 
ylim([.9995*min(Ypol) 1.0001*max(Ypol)])


end
% 
% figure(3)
% plot2Dgrid(lac_x,lac_y)
% title 'Curvilinear Grid (Aerial View)'
% 
% 
% figure(3)
% plot2Dgrid(ihs_hor,ihs_vert)
% title 'Curvilinear Grid (Sectional View)'
% hold off




% DISPLAY 3D GRID
figure(4)
s = surf(myfinalx,myfinaly,myfinalz);colormap(jet); hold off ;
s.EdgeColor = 'k';
pbaspect([1.5 1 .25 ])
lim_min_y = min(Resy); %min_x = min(lim_minx); 
lim_max_y = max(Resy); %lim_max_x = max(lim_max_x);
if CONCAVE_DIR ==1
a1 = 1.0001*lim_min_y(2); a2 = lim_max_y(1);
ylim([a1 a2]); zlim([0 PBH])
elseif CONCAVE_DIR ==2
a1 = lim_min_y(2); a2 = 0.9999*lim_max_y(1);
ylim([a1 a2]); zlim([0 PBH])
end
title '3D Point Bar Curvilinear Grid'


figure
surf(reshape(xx_beg,blks_alongcurve+1,blks_along_sigmoid+1),reshape(yy_beg,blks_alongcurve+1,blks_along_sigmoid+1),reshape(zz_beg,blks_alongcurve+1,blks_along_sigmoid+1));
hold on; colormap(jet);pbaspect([aspx aspy aspz ]) %shading interp;  
%s.EdgeColor = 'none';
xlabel('x-coordinates'); ylabel('y-coordinates'); zlabel('z-coordinates')
surf(reshape(xx_end,blks_alongcurve+1,blks_along_sigmoid+1),reshape(yy_end,blks_alongcurve+1,blks_along_sigmoid+1),reshape(zz_end,blks_alongcurve+1,blks_along_sigmoid+1))
hold on;colormap(jet);pbaspect([aspx aspy aspz ]) %shading interp;
xlabel('x-coordinates'); ylabel('y-coordinates'); zlabel('z-coordinates')

surf(reshape(XX_beg,N_B+1,blks_along_sigmoid+1),reshape(YY_beg,N_B+1,blks_along_sigmoid+1),reshape(ZZ_end,N_B+1,blks_along_sigmoid+1))
hold on; colormap(jet); pbaspect([aspx aspy aspz ]) %shading interp;
xlabel('x-coordinates'); ylabel('y-coordinates'); zlabel('z-coordinates')

surf(reshape(XX_end,N_B+1,blks_along_sigmoid+1),reshape(YY_end,N_B+1,blks_along_sigmoid+1),reshape(ZZ_end,N_B+1,blks_along_sigmoid+1))
xlabel('x-coordinates'); ylabel('y-coordinates'); zlabel('z-coordinates')
hold on; colormap(jet); pbaspect([aspx aspy aspz ]) % shading interp;
surf(reshape(xx_down,blks_alongcurve+1,N_B+1),reshape(yy_down,blks_alongcurve+1,N_B+1),reshape(zz_down,blks_alongcurve+1,N_B+1))
xlabel('x-coordinates'); ylabel('y-coordinates'); zlabel('z-coordinates')
hold on; colormap(jet); pbaspect([aspx aspy aspz ]) %shading interp;
surf(reshape(xx_up,blks_alongcurve+1,N_B+1),reshape(yy_up,blks_alongcurve+1,N_B+1),reshape(zz_up,blks_alongcurve+1,N_B+1))
colormap(jet);pbaspect([aspx aspy aspz ]) %shading interp;
xlabel('x-coordinates'); ylabel('y-coordinates'); zlabel('Point Bar Height')
title 'Equivalent Rectilinear Grid'
pbaspect([1 1 .2])
hold off







end