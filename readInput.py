from typing import TextIO, Tuple, Union
import numpy as np

import headLine
import readMatrix


def readInput(
    iread: TextIO,
    id_char: Union[str, int],
    nnod: int,
    nco: int,
    ndn: int,
    nbc: int,
    nmat: int,
    nsec: int,
    iforce: int,
    nde: int,
    itp: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Ported from MATLAB: INPUT.m

    Reads nodal coordinates, boundary conditions, external nodal loads,
    element connectivity (IDBC), local y-axis direction cosines (VECTY, space frame only),
    fixed-end forces (FEF, when IFORCE=2), material properties (PROP), and
    cross-sectional properties (SECT).

    Parameters
    ----------
    iread : TextIO
        Open text file handle to the input source.
    id_char : str | int
        Identifier character that starts each data block line in the input file.
        Behavior matches HEADLINE.m: we scan forward until a line whose first char equals `id_char`.
    nnod : int
        Number of nodes.
    nco : int
        Number of coordinates per node.
    ndn : int
        Number of DOFs per node.
    nbc : int
        Number of beam-column elements.
    nmat : int
        Number of material types.
    nsec : int
        Number of cross-sectional types.
    iforce : int
        1 -> only concentrated nodal loads (no FEF block).
        2 -> need to read Fixed-End Forces (FEF) for each member.
    nde : int
        Number of DOFs per element.
    itp : int
        Frame type. If 6 (space frame), VECTY block is expected and read.

    Returns
    -------
    COOR, NFIX, EXLD, FEF, IDBC, PROP, SECT, VECTY : tuple of np.ndarray
        - COOR shape (nco, nnod), float
        - NFIX shape (ndn, nnod), int  (values -1, 0, or node-number for double-node technique)
        - EXLD shape (ndn, nnod), float
        - FEF  shape (nde, nbc), float; empty (0 * 0) when IFORCE==1
        - IDBC shape (5, nbc),  int  ([1]=node1, [2]=node2, [3]=mat id, [4]=sect id, [5]=omitted)
        - PROP shape (5, nmat), float
        - SECT shape (5, nsec), float
        - VECTY shape (3, nbc), float; empty (0 * 0) when ITP!=6

    Notes
    -----
    - Parsing uses `ReadMatrix.ReadMatrix`, which for each matrix block reads `col` lines and
      takes tokens 2..(row+1) per line (skipping the 1st token), mirroring MATLAB `num(2:row+1)`.
    - The printed echoes mimic the MATLAB script (`disp(...)`).
    - IDBC and NFIX are returned as integer arrays; other blocks are float arrays.
    - Node numbering inside IDBC is preserved as found in file (typically 1-based).
    """
    # COOR - Nodal coordinates
    headLine.headline(id_char, iread)
    COOR = readMatrix.readMatrix(iread, nco, nnod)
    print("COOR :")
    print(COOR.T)

    # NFIX - DOF flags per node
    headLine.headline(id_char, iread)
    NFIX = readMatrix.readMatrix(iread, ndn, nnod).astype(int, copy=False)
    print("NFIX :")
    print(NFIX.T)

    # EXLD - external nodal loads
    headLine.headline(id_char, iread)
    EXLD = readMatrix.readMatrix(iread, ndn, nnod)
    print("EXLD :")
    print(EXLD.T)

    # IDBC - element identification (connectivity & types)
    headLine.headline(id_char, iread)
    IDBC = readMatrix.readMatrix(iread, 5, nbc).astype(int, copy=False)
    print("IDBC :")
    print(IDBC.T)

    # VECTY - local y-axis direction cosines (space frame only)
    if itp != 6:
        VECTY = np.empty((0, 0), dtype=float)
    else:
        headLine.headline(id_char, iread)
        VECTY = readMatrix.readMatrix(iread, 3, nbc)
        print("VECTY:")
        print(VECTY.T)

    # FEF - fixed-end forces (when IFORCE==2)
    if iforce == 1:
        FEF = np.empty((0, 0), dtype=float)
    else:
        headLine.headline(id_char, iread)
        FEF = readMatrix.readMatrix(iread, nde, nbc)
        print("FEF:")
        print(FEF.T)

    # PROP - material properties
    headLine.headline(id_char, iread)
    PROP = readMatrix.readMatrix(iread, 5, nmat)
    print("PROP:")
    print(PROP.T)

    # SECT - section properties
    headLine.headline(id_char, iread)
    SECT = readMatrix.readMatrix(iread, 5, nsec)
    print("SECT:")
    print(SECT.T)

    return COOR, NFIX, EXLD, FEF, IDBC, PROP, SECT, VECTY