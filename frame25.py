from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import headLine
import readInput
import drawStructure


@dataclass
class Frame25Result:
    title: str
    funit: str
    lunit: str
    itp: int
    nco: int
    ndn: int
    nne: int
    nde: int
    nnod: int
    nbc: int
    nmat: int
    nsec: int
    coor: np.ndarray
    nfix: np.ndarray
    exld: np.ndarray
    fef: np.ndarray
    idbc: np.ndarray
    prop: np.ndarray
    sect: np.ndarray
    vecty: np.ndarray
    figure_path: Optional[Path]


FTYPE = ("BEAM", "PLANE TRUSS", "PLANE FRAME", "PLANE GRID", "SPACE TRUSS", "SPACE FRAME")
IPR = np.array([[1, 2, 2, 2, 3, 3],
                [2, 2, 3, 3, 3, 6]], dtype=int)


def _parse_ints_from_line(line: str, n: int) -> Tuple[int, ...]:
    toks = line.strip().replace(",", " ").split()
    if len(toks) < n:
        raise ValueError(f"Expected at least {n} numbers, got {len(toks)}: {line!r}")
    vals = []
    for k in range(n):
        x = float(toks[k])
        vals.append(int(round(x)))
    return tuple(vals)


def frame25(filename: str | Path, save_figure: bool = True) -> Frame25Result:

    start_time = datetime.now()

    filename = Path(filename)
    ipt_path = filename.with_suffix(".ipt")

    if not ipt_path.exists():
        raise FileNotFoundError(f"Input file not found: {ipt_path}")

    with ipt_path.open("r", encoding="utf-8") as f:
        title = (f.readline() or "").rstrip("\n")
        funit = (f.readline() or "").strip()
        lunit = (f.readline() or "").strip()

        ID = "*"
        headLine.headline(ID, f)
        args_line = f.readline()
        NNOD, NBC, NMAT, NSEC, ITP, NNE, IFORCE = _parse_ints_from_line(args_line, 7)

        if not (1 <= ITP <= 6):
            raise ValueError(f"ITP must be 1..6, got {ITP}")
        nco = int(IPR[0, ITP - 1])
        ndn = int(IPR[1, ITP - 1])
        nde = ndn * NNE

        COOR, NFIX, EXLD, FEF, IDBC, PROP, SECT, VECTY = readInput.readInput(
            iread=f,
            id_char=ID,
            nnod=NNOD,
            nco=nco,
            ndn=ndn,
            nbc=NBC,
            nmat=NMAT,
            nsec=NSEC,
            iforce=IFORCE,
            nde=nde,
            itp=ITP,
        )

    fmt = "k"
    fig, ax = drawStructure.drawStructure(ITP, COOR, IDBC, NBC, lunit, fmt)

    fig.tight_layout()
    out_pdf: Optional[Path] = None
    if save_figure:
        out_pdf = filename.with_name(f"{filename}-draw.pdf")
        fig.savefig(out_pdf, format="pdf", dpi=300)

    plt.close(fig)

    end_time = datetime.now()
    duration = end_time - start_time

    print(f"Start Time: {start_time}")
    print(f"End   Time: {end_time}")
    print(f"Duration  : {duration}")

    return Frame25Result(
        title=title,
        funit=funit,
        lunit=lunit,
        itp=ITP,
        nco=nco,
        ndn=ndn,
        nne=NNE,
        nde=nde,
        nnod=NNOD,
        nbc=NBC,
        nmat=NMAT,
        nsec=NSEC,
        coor=COOR,
        nfix=NFIX,
        exld=EXLD,
        fef=FEF,
        idbc=IDBC,
        prop=PROP,
        sect=SECT,
        vecty=VECTY,
        figure_path=out_pdf,
    )