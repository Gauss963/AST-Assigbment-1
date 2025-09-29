function FRAME25(FILENAME)
%FRAME25 - A linear analysis program for framed structures
%
% Syntax:
%   FRAME25(FILENAME)
%
% Description:
%   Reads {FILENAME}.ipt as input for a framed-structure linear analysis.
%   The program parses nodes/elements/materials/sections, handles loading
%   options (concentrated vs. fixed-end forces), and can draw the geometry.
%
% Programmer:
%   Gauss Chang (R14521220)
%   Supervised by Assistant Professor Chun-Yu Ke
%   For the course: Advanced Structural Theory
%   Department of Civil Engineering, National Taiwan University
%   Term: Fall 2025  © All Rights Reserved
%
% Inputs:
%   FILENAME  - (char|string) Base name (without extension) of the input
%               file. The program expects:
%                 - {FILENAME}.ipt : input data file
%                 - {FILENAME}.dat : output listing (echoed input, nodal
%                                    displacements, member end forces)
%
% Variables (internal):
%   NNOD   = number of nodes
%   NBC    = number of Beam-column elements
%   NCO    = number of coordinates per node
%   NDN    = number of DOFs per node
%   NNE    = number of nodes per element
%   NDE    = number of DOFs per element
%   NMAT   = number of material types
%   NSEC   = number of cross-sectional types
%   IFORCE = loading flag:
%              1 -> only concentrated loads are applied
%              2 -> fixed-end forces are required
%                   (e.g., distributed loads, fabrication errors,
%                    temperature change, support settlement)
%
% Strings from input (legacy note "CHARACTERS"):
%   FUNIT  = unit of force (e.g., kN, kip)
%   LUNIT  = unit of length (e.g., mm, in)
%
% Frame type legend:
%   FRAME TYPE      ITP   NCO   NDN   (NCO and NDN are stored in array IPR)
%   BEAM             1     1     2
%   PLANAR TRUSS     2     2     2
%   PLANAR FRAME     3     2     3
%   PLANAR GRID      4     2     3
%   SPACE TRUSS      5     3     3
%   SPACE FRAME      6     3     6
%
% Notes:
%   - Place this help block immediately after the function line (no blank
%     lines in between). The first line "%FRAME25 - ..." is the H1 line
%     used by HELP/LOOKFOR/Doc and by VS Code hover.
%   - Each subfunction in this file may include its own help block using
%     the same style, directly below its own function line.
%
% See also:
%   INPUT, drawStructure, exportgraphics

FTYPE = {'BEAM'; 'PLANE TRUSS'; 'PLANE FRAME'; 'PLANE GRID'; 'SPACE TRUSS'; 'SPACE FRAME'};
IPR = [1, 2, 2, 2, 3, 3; 2, 2, 3, 3, 3, 6];

if nargin == 0
    [~,FILENAME] = fileparts(uigetfile('*.ipt','Select Input file'));
end

% Get starting time
startTime = datetime("now");

% {FILENAME}.ipt is the input data file.
% {FILENAME}.dat includes the output of the input data and the nodal 
% displacements and member end forces.

IREAD = fopen([FILENAME '.ipt'], 'r');

% Read in the problem title and the structural data
TITLE = fgets(IREAD);
FUNIT = fgets(IREAD);
LUNIT = fgets(IREAD);
ID = '*';
HEADLINE(ID,IREAD);
line = fgets(IREAD);
args = str2num(line);
temp = num2cell(args(1:7));
[NNOD, NBC, NMAT, NSEC, ITP, NNE, IFORCE] = deal(temp{:});
NCO = IPR(1, ITP);
NDN = IPR(2, ITP);
NDE = NDN * NNE;

% Read the remaining data
[COOR, NFIX, EXLD, FEF, IDBC, PROP, SECT, VECTY] = INPUT(IREAD, ID, NNOD, NCO, NDN, NBC, NMAT, NSEC, IFORCE,NDE, ITP);
FORMAT='k';
drawStructure(ITP,COOR,IDBC,NBC,LUNIT,FORMAT);
% exportgraphics(gcf, 'Problem-1-structure.pdf', 'ContentType', 'vector');
fig = gcf; ax = gca;
set(fig,'Color','w'); 
set(ax,'Color','w', 'XColor','k','YColor','k','ZColor','k');
set(findall(fig,'Type','text'),'Color','k');
exportgraphics(fig,'structure.pdf', 'ContentType','vector', 'BackgroundColor','white');
fclose(IREAD);
endTime = datetime("now");

duration = endTime - startTime;
disp(['Start Time: ', char(startTime)]);
disp(['End   Time: ', char(endTime)]);
disp(['Duration   : ', char(duration)]);

end