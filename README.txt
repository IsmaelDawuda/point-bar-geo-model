

This source code  recreates the channel meander path ( conditioned to well data) and possible point bar well locations 
to  model a gridded point bar geometry, and assign its petrophysical property distributions obtained from Geostatistical simulation.

The coordinates of wells should be in feets. 

Download all data and run the 'point-bar-geomodel.m' funtion using the "Demo.m" file

 ---- RUNNIG THE CODE FROM THE DEMO FILE ---

The code in this file activates the 'point-bar-geomodel.m' function
To run the demo file, the following input parameters  must be specified. 
filename.txt    :       % the file containing the coordinates of the channel locations; the first and second columns are the x and y coordinates, respectively   

load loc_wells.txt;     % the file containing the coordinates of wells (like injection wells, producers, observation wells etc ) that may be present in the reservoir;  the first and second columns are the x and y coordinates, respectively
                        % Specifying this ensures that the relative positions  of the wells in the reservoir are preserved during the channel migration or geometry peturbation.  
                        % if there is no well present, leave the content of this file, it will not affect the final results because the coordinates of the point bar grid nodes, and the location of these wells will be reported separately;  
                        % Please note that you may change the contents of these files but the file names must be retained

nx             :        % number of grid blocks along the x-axis 
ny             :        % number of grid blocks along the y-axis 
nz             :        % number of grid blocks along the z-axis 
PBH            :        % point bar height in ft
DTOP           :        % depth to top of formation, in ft

property_SGSIM :        %the property data which has been modeled using GEOSTATISTICS (SGSIM and ELLIPSIM )

Once the above input parameters are specified, the 'Demo' file can be run

+++++++The following prompt messages will appear once the code begins to run+++++

-----RESPONSE TO PROMPT MESSAGES ---
A) Enter a rotation angle to align points with the global X-Y plane, Enter 0 if points are already aligned:
---The Code assumes that the channel is progressing in the East-West Direction; this message prompts you to enter 
a rotation angle to align points in the E-W direction. At the  end of the computations, the points will be back transformed
to their original orientation.  For this particular example, the coordinates are already aligned, so you may enter 0.


B) Which of the  concave parts should be used? the number should be <= number of concave parts: 
---The Meander paths may have more than one bends (i.e., more than one concave regions) depending on the data input; 
this message allows you specify which of the concave parts to use to model the point bar
For this particular dataset example,  there are three concave parts; so you may enter either 1, 2 or 3


C) What is the angular displacement of this channel? :
The Channel path is approximated with a Sine Generation Function as proposed by Langbein and Leopold, 1966.
This requires that an angular displacment is specified to define a unique channel path according to the input data provided.
For this particular dataset example,  you may enter any angle within a range of 25 to 35. 

Upon responding to the prompt messages, the 'point-bar-geo-model.m' code would be executed, and 
the final gridded point bar geometry will be saved  as a text file in an ascii format  as 'corner_point_grid.txt'
if there are any well specified, the locations will also be be reported and saved as a textile in an ascii format as 'well_locs_locations.xt'

In addition, the follwoing would be displayed graphically:

1) the raw data showing the channel meander locations
2) channel meander locations sorted in the direction of channel progression
3) Randomly generated possible point bar well locations and  Channel path conditioned to go through channel nodes and bending to accomodate point bar locations on the concave side of the meander
1) Grid for the Lateral Accretions in 2D
2) Grid for the inclined heterolithic stratifications  in 2D
3) 3D curvilinear and rectilinear Grid for the entire point Bar
4) Property distribution for the Point Bar in 3D




Please note:
For better visualization, you may consider exporting the grid and the reservoir properties to a different software with better graphic quality like CMG, PETREL etc.
