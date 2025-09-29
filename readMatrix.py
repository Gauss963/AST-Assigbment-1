from typing import TextIO
import numpy as np


def readMatrix(iread: TextIO, row: int, col: int) -> np.ndarray:
    mat = np.zeros((row, col), dtype = float)

    for j in range(col):
        line = iread.readline()
        if line == "":
            raise EOFError(f"Unexpected end of file while reading column {j+1}/{col}.")

        tokens = line.strip().replace(",", " ").split()

        if len(tokens) < row + 1:
            raise ValueError(
                f"Line {j+1} needs at least {row+1} numbers, got {len(tokens)}: {line!r}"
            )

        try:
            values = [float(tok) for tok in tokens[1 : row + 1]]
        except ValueError as e:
            raise ValueError(f"Non-numeric token on line {j+1}: {e}")

        mat[:, j] = values

    return mat