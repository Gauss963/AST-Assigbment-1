function [COOR,NFIX,EXLD,FEF,IDBC,PROP,SECT,VECTY] = INPUT(IREAD,ID,NNOD,NCO,NDN,NBC,NMAT,NSEC,IFORCE,NDE,ITP)

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

% COOR - Nodal coordinates
HEADLINE(ID,IREAD);
COOR = ReadMatrix(IREAD,NCO,NNOD);
disp('COOR :');
disp(COOR');

% NFIX - DOF in each nodes
HEADLINE(ID,IREAD);
NFIX = ReadMatrix(IREAD,NDN,NNOD);
disp('NFIX :');
disp(NFIX');

% EXLD - external load
HEADLINE(ID,IREAD);
EXLD = ReadMatrix(IREAD,NDN,NNOD);
disp('EXLD :');
disp(EXLD');

% IDBC - Nodes in each element
HEADLINE(ID,IREAD);
IDBC = ReadMatrix(IREAD,5,NBC);
disp('IDBC :');
disp(IDBC');

if ITP ~=6
    VECTY = [];
else
    HEADLINE(ID,IREAD);
    VECTY = ReadMatrix(IREAD,3,NBC);
    disp('VECTY: ');
    disp(VECTY');
end

if IFORCE ==1
    FEF = [];
else
    HEADLINE(ID,IREAD);
    FEF = ReadMatrix(IREAD,NDE,NBC);
    disp('FEF: ');
    disp(FEF');
end

HEADLINE(ID,IREAD);
PROP = ReadMatrix(IREAD,5,NMAT);
disp('PROP: ');
disp(PROP');

HEADLINE(ID,IREAD);
SECT = ReadMatrix(IREAD,5,NSEC);
disp('SECT: ');
disp(SECT');


end
