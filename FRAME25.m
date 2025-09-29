%     Advanced Structural Theory 2025
%
%     Program Assignment No. 1 (Weight: 1)
%
%     Note: Each program assignment has a specific weight. Generally, the 
%     higher the weight, the more time you are expected to dedicate to the 
%     assignment.
%
%        Assigned:  09/25/2025
%        Due:       10/09/2025
%
%      (1) Complete the function INPUT.
%      (2) Implement the main program FRAME25 up to the point marked
%          % ^^ Up to this point --- PROG 1 ^^
%      (3) Test the program using the provided problem in "prog1.pdf". You
%          will need to create an input file and run the program, ensuring
%          that the output data (from function INPUT) matches the input 
%          data. Additionally, use the drawStructure function to verify 
%          the geometry of the structures.
%      (4) Submit the following in a compressed archive (*.zip) to COOL 
%          with your student ID as the file name:
%          (a) Program source code "FRAME25.m"
%          (b) Input files "*.ipt"
%
function FRAME25(FILENAME)
% FRAME25: A linear analysis program for framed structures
%..........................................................................
%    Programmer:  YOUR NAMES and STUDENT IDs(both your and your partner's)
%                 Supervised by Assistant Professor Chun-Yu Ke
%                 For the course: Advanced Structural Theory
%                 Department of Civil Engineering
%                 National Taiwan University
%                 Fall 2025 © All Rights Reserved
%..........................................................................
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
%..........................................................................
%    FRAME TYPE    ITP  NCO  NDN   (NCO and NDN are stored in Array IPR)
%    BEAM           1    1    2
%    PLANAR TRUSS   2    2    2
%    PLANAR FRAME   3    2    3
%    PLANAR GRID    4    2    3
%    SPACE  TRUSS   5    3    3
%    SPACE  FRAME   6    3    6

FTYPE = {'BEAM';'PLANE TRUSS';'PLANE FRAME';'PLANE GRID';'SPACE TRUSS';'SPACE FRAME'};
IPR = [1,2,2,2,3,3;2,2,3,3,3,6];

if nargin == 0 % no input argument
    % Open file with user interface
    [~,FILENAME] = fileparts(uigetfile('*.ipt','Select Input file'));
end

% Get starting time
startTime = clock;

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
[NNOD,NBC,NMAT,NSEC,ITP,NNE,IFORCE] = deal(temp{:});
NCO = IPR(1,ITP);
NDN = IPR(2,ITP);
NDE = NDN*NNE;

% Read the remaining data
%%%%% Finish the INPUT function



% ˇ ˇ Starting from here --- PROG 1

% [COOR,~] = INPUT(IREAD,ID,NNOD,NCO,~); % you need to add variables here
% please use the function - drawStructure - here to plot the structure
fclose(IREAD);

% ^ ^ Up to Here --- PROG 1



% 
% % DOF numbering
% % Finish the IDMAT function.
% [IDND,NEQ] = IDMAT( ~ );
% 
% % Create the member DOF table:  LM(NDE,NBC)
% % Finish the MEMDOF function.
% LM = MEMDOF( ~ );
% 
% % Compute the semi-band width,NSBAND, of the global stiffness matrix.
% % Finish the SEMIBAND function.
% NSBAND = SEMIBAND( ~ );
% 
% % Form the global load vector GLOAD(NEQ) from the concentrated nodal loads.
% % Finish the LOAD function.
% GLOAD = LOAD(EXLD, ~ );
% 
% % ^^ UP TO HERE --- PROG 2 ^^
% 
% % Form the global stiffness matrix GLK(NEQ,NSBAND) and obtain the
% % equivalent nodal vector by assembling -(fixed-end forces) of each member
% % into the load vector.
% [GLK,GLOAD] = FORMKP(COOR,IDBC,VECTY,PROP,SECT,LM,FEF,GLOAD,NNOD,NBC,NMAT...
%     ,NSEC,IFORCE,ITP,NCO,NDN,NDE,NNE,NEQ);
% 
% % ^^ UP TO HERE --- PROG 3 ^^
% 
% DISP = SOLVE(GLK, ~ );
% 
% % Determine the member end forces ELFOR(NDE,NBC)
% ELFOR = FORCE( ~ );
% 
% % Get ending time and count the elapased time
% endTime = clock;
% 
% % Print out the results
% IWRITE = fopen([FILENAME '.dat'], 'w');
% OUTPUT( ~ );
% fclose(IWRITE);
% 
% 
% % ^^ UP TO HERE --- PROG 4 ^^

end

function HEADLINE(ID,IREAD)
    while ~feof(IREAD)
        temp = fgets(IREAD);
        if ~isempty(temp) && temp(1) == ID
            return
        end
    end
end

function [COOR ...
   ] = INPUT(IREAD,ID,NNOD,NCO ...
   )

%..........................................................................
%
%    PURPOSE: Reads in the following data: nodal coordinates,
%             boundary conditions, external nodal loads, identification
%             data (connectivity), direction cosines of the local y-axis
%             (for space frame only; ITP=6), fixed-end forces, material
%             and cross sectional properties.
%
%    VARIABLES:
%      IREAD          = index of input stream
%      ID             = identifier charactor of input file
%      NNOD           = number of nodes
%      NBC            = number of Beam-column elements
%      NMAT           = number of material types
%      NSEC           = number of cross-sectional types
%      ITP            = frame tyoe
%      NCO            = number of coordinates per node
%      NDN            = number of DOFs per node
%      NDE            = number of DOFs per element
%      IFORCE         = indicator for fixed-end forces (FEF)
%                     = 1 for the case of ONLY concentrated nodal loads;
%                       no need to input FEF.
%                     = 2 for the following cases: distributed loads,
%                       temperature change, support settlement;
%                       need to input FEF for each member.
%      COOR(NCO,NNOD) = nodal coordinates
%      NFIX(NDN,NNOD) = Boolian id matrix for specifying boundary conditions
%                     = -1 for restrained d.o.f.
%                     =  0 for free d.o.f.
%                     = the first node number for the d.o.f.
%                       that has the same d.o.f. of the second node
%                       when "the double node technique" is used.
%
%      EXLD(NDN,NNOD) = external nodal loads
%      IDBC(5,NBC)    = Beam-column identification data
%                       (1,*) = local node 1
%                       (2,*) = local node 2
%                       (3,*) = material type.
%                       (4,*) = section type.
%                       (5,*) = omitted.
%      VECTY(3,NBC)   = direction cosines of the local y-axis
%                       VECTY is needed only when ITP=6 (space frame)
%      FEF(NDE,NBC)   = FEF is needed only when IFORCE=2
%      PROP(5,NMAT)   = material properties
%                       (1,*) = Young's modulus
%                       (2,*) = Poision ratio.
%                       (3,*) = omitted.
%                       (4,*) = omitted.
%                       (5,*) = omitted.
%      SECT(5,NSEC)   = Beam-column properties
%                       (1,*) = cross-sectional area A.
%                       (2,*) = moment of inertia Iz (3D)
%                       (3,*) = moment of inertia Iy (3D)
%                       Note that for 2D problems, only Iz is required.
%                       (4,*) = torsional constant J.
%                       (5,*) = omitted.
%..........................................................................




% ˇ ˇ Starting from here --- PROG 1

% COOR - Nodal coordinates
HEADLINE(ID,IREAD);
COOR = ReadMatrix(IREAD,NCO,NNOD);
disp('COOR :');
disp(COOR');
% Please finish the rest of the function here.
% You have to display all the matrix in .ipt to the command window.(COOR,NFIX,EXLD,IDBC,PROP,SECT)
% ...???
% ...???


% ^ ^ Up to here --- PROG 1
end






function mat = ReadMatrix(IREAD, row, col)
% for INPUT function
    mat = zeros(row,col);
    for j = 1:col
        line = fgets(IREAD);
        num = str2num(line);
        mat(:,j) = num(2:row+1);
    end
end

function drawStructure(ITP,COOR,IDBC,NBC,LUNIT,FORMAT)
% draw the structure
    switch ITP
        case 1
            for e=1:NBC
                plot([COOR(IDBC(1,e)),COOR(IDBC(2,e))],[0,0],FORMAT,'linewidth',2)
                xlabel(['X ', LUNIT])
                title('Beam')
                hold on
            end
            hold off
        case {2,3}
            for e=1:NBC
                plot([COOR(1,IDBC(1,e)),COOR(1,IDBC(2,e))],...
                    [COOR(2,IDBC(1,e)),COOR(2,IDBC(2,e))],FORMAT,'linewidth',2)
                xlabel(['X ', LUNIT])
                ylabel(['Y ', LUNIT])
                if ITP==2
                    title('2D truss')
                else
                    title('2D frame')
                end
                hold on
            end
            axis equal
            hold off
        case 4
            for e=1:NBC
                plot([COOR(1,IDBC(1,e)),COOR(1,IDBC(2,e))],...
                    [COOR(2,IDBC(1,e)),COOR(2,IDBC(2,e))],FORMAT,'linewidth',2)
                xlabel(['X ', LUNIT])
                ylabel(['Z ', LUNIT])
                title('grid')
                hold on
            end
            axis equal
            hold off
        case {5,6}
            for e=1:NBC
                plot3([COOR(1,IDBC(1,e)),COOR(1,IDBC(2,e))],...
                    [COOR(2,IDBC(1,e)),COOR(2,IDBC(2,e))],...
                    [COOR(3,IDBC(1,e)),COOR(3,IDBC(2,e))],FORMAT,'linewidth',2)
                xlabel(['X ', LUNIT])
                ylabel(['Y ', LUNIT])
                zlabel(['Z ', LUNIT])
                if ITP==5
                    title('3D truss')
                else
                    title('3D frame')
                end
                hold on
            end
            axis equal
            hold off
    end
end