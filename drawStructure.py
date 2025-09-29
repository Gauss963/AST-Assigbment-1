from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def _set_3d_axes_equal(ax) -> None:
    """Make 3D axes have equal scale—like MATLAB's axis equal for 3D."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]
    y_range = y_limits[1] - y_limits[0]
    z_range = z_limits[1] - z_limits[0]
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    x_mid = 0.5 * (x_limits[0] + x_limits[1])
    y_mid = 0.5 * (y_limits[0] + y_limits[1])
    z_mid = 0.5 * (z_limits[0] + z_limits[1])

    ax.set_xlim3d([x_mid - plot_radius, x_mid + plot_radius])
    ax.set_ylim3d([y_mid - plot_radius, y_mid + plot_radius])
    ax.set_zlim3d([z_mid - plot_radius, z_mid + plot_radius])


def drawStructure(
    itp: int,
    coor: np.ndarray,
    idbc: np.ndarray,
    nbc: int,
    lunit: str,
    fmt: str,
    ax: Optional[plt.Axes] = None,
    linewidth: float = 2.0,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Ported from MATLAB: drawStructure.m

    Parameters
    ----------
    itp : int
        結構型別：1=beam(1D), 2=2D truss, 3=2D frame, 4=grid (X-Z),
                  5=3D truss, 6=3D frame。
    coor : np.ndarray
        節點座標，形狀預期為 (NCO, NNOD)。例如 2D: (2, NNOD)、3D: (3, NNOD)。
    idbc : np.ndarray
        元素識別資料，形狀 (5, NBC)。其中第 1、2 列為兩端節點編號（通常為 1-based）。
    nbc : int
        元素數（loop 次數）。若與 idbc.shape[1] 不同，仍以此參數為主。
    lunit : str
        座標單位字串（e.g., 'm', 'cm'）。
    fmt : str
        matplotlib 線條/標記格式（e.g. '-k', '--o'）。
    ax : Optional[plt.Axes]
        既有的 Axes。若為 None 則自動建立。
    linewidth : float
        線寬，預設 2.0。

    Returns
    -------
    (fig, ax) : Tuple[Figure, Axes]
        供呼叫端後續保存或調整。

    Notes
    -----
    - 依 MATLAB 寫法，節點編號採用檔內的 1-based；本函式在索引時會自動轉成 0-based。
    - 2D 模式使用 `ax.set_aspect('equal', adjustable='box')`；3D 模式用自訂 `_set_3d_axes_equal`。
    """
    idbc = np.asarray(idbc)
    coor = np.asarray(coor)

    # 建立 axes
    is3d = itp in (5, 6)
    if ax is None:
        fig = plt.figure()
        if is3d:
            ax = fig.add_subplot(111, projection="3d")
        else:
            ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    # 取方便的 coor 索引（coor 形狀為 (NCO, NNOD)）
    def _x(idx0: int) -> float:
        return float(coor[0, idx0])

    def _y(idx0: int) -> float:
        return float(coor[1, idx0])

    def _z(idx0: int) -> float:
        return float(coor[2, idx0])

    # 逐元素畫線
    for e in range(nbc):
        # 讀取兩端節點（1-based -> 0-based）
        i = int(idbc[0, e]) - 1
        j = int(idbc[1, e]) - 1

        if itp == 1:
            # Beam: X 軸畫在 y=0
            xs = [_x(i), _x(j)]
            ys = [0.0, 0.0]
            ax.plot(xs, ys, fmt, linewidth=linewidth)

        elif itp in (2, 3):
            # 2D: X-Y
            xs = [_x(i), _x(j)]
            ys = [_y(i), _y(j)]
            ax.plot(xs, ys, fmt, linewidth=linewidth)

        elif itp == 4:
            # grid: X-Z
            xs = [_x(i), _x(j)]
            zs = [_z(i), _z(j)]
            ax.plot(xs, zs, fmt, linewidth=linewidth)
        elif itp in (5, 6):
            # 3D: X-Y-Z
            xs = [_x(i), _x(j)]
            ys = [_y(i), _y(j)]
            zs = [_z(i), _z(j)]
            ax.plot(xs, ys, zs, fmt, linewidth=linewidth)
        else:
            raise ValueError(f"Unknown ITP={itp}")

    # 標籤與標題
    if itp == 1:
        ax.set_xlabel(f"X {lunit}")
        ax.set_title("Beam")
    elif itp in (2, 3):
        ax.set_xlabel(f"X {lunit}")
        ax.set_ylabel(f"Y {lunit}")
        ax.set_title("2D truss" if itp == 2 else "2D frame")
        ax.set_aspect("equal", adjustable="box")
    elif itp == 4:
        ax.set_xlabel(f"X {lunit}")
        ax.set_ylabel(f"Z {lunit}")
        ax.set_title("grid")
        ax.set_aspect("equal", adjustable="box")
    elif itp in (5, 6):
        ax.set_xlabel(f"X {lunit}")
        ax.set_ylabel(f"Y {lunit}")
        ax.set_zlabel(f"Z {lunit}")
        ax.set_title("3D truss" if itp == 5 else "3D frame")
        _set_3d_axes_equal(ax)

    return fig, ax