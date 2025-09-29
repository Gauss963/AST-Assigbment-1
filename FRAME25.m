function FRAME25(FILENAME)
% FRAME25: A linear analysis program for framed structures
%--------------------------------------------------------------------------
%    Programmer:  Gauss Chang (R14521220)
%                 Supervised by Assistant Professor Chun-Yu Ke
%                 For the course: Advanced Structural Theory
%                 Department of Civil Engineering
%                 National Taiwan University
%                 Fall 2025 © All Rights Reserved
%--------------------------------------------------------------------------
%    VARIABLES:
%        NNOD   = number of nodes
%        NBC    = number of Beam-column elements
%        NCO    = number of coordinates per node
%        NDN    = number of DOFs per node
%        NNE    = number of nodes per element
%        NDE    = number of DOFs per element
%        NMAT   = number of material types
%        NSEC   = number of cross-sectional types
%        IFORCE = 1 if only concentrated loads are applied
%               = 2 if fixed-end forces are required.
%                   (e.g. problems with distributed loads, fabrication
%                   errors, temperature change, or support settlement)
%    CHARACTERS
%        FUNIT  = unit of force (such as kN and kip)
%        LUNIT  = unit of length (such as mm and in)
%--------------------------------------------------------------------------
%    FRAME TYPE    ITP  NCO  NDN   (NCO and NDN are stored in Array IPR)
%    BEAM           1    1    2
%    PLANAR TRUSS   2    2    2
%    PLANAR FRAME   3    2    3
%    PLANAR GRID    4    2    3
%    SPACE  TRUSS   5    3    3
%    SPACE  FRAME   6    3    6

FTYPE = {'BEAM';'PLANE TRUSS';'PLANE FRAME';'PLANE GRID';'SPACE TRUSS';'SPACE FRAME'};
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

end