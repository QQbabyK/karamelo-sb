import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_log(file_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """读取 log.mpm，并返回时间和悬臂梁端部位移。"""
    if not file_path.exists():
        raise FileNotFoundError(f"找不到结果文件：{file_path}")

    data = np.genfromtxt(
        file_path,
        names=True,
        dtype=float,
        encoding="utf-8",
    )

    required_columns = {"Time", "dy_tip"}
    if data.dtype.names is None:
        raise ValueError("无法识别 log.mpm 的表头。")

    missing = required_columns.difference(data.dtype.names)
    if missing:
        raise ValueError(
            f"log.mpm 缺少列：{', '.join(sorted(missing))}；"
            f"现有列：{', '.join(data.dtype.names)}"
        )

    time = np.atleast_1d(data["Time"])
    displacement = np.atleast_1d(data["dy_tip"])

    valid = np.isfinite(time) & np.isfinite(displacement)
    time = time[valid]
    displacement = displacement[valid]

    if time.size == 0:
        raise ValueError("log.mpm 中没有有效数据。")

    return time, displacement


def main() -> None:
    parser = argparse.ArgumentParser(
        description="绘制 Karamelo 悬臂梁自由端位移随时间变化曲线。"
    )
    parser.add_argument(
        "-i",
        "--input",
        default="log.mpm",
        help="输入日志文件，默认：log.mpm",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="dy_tip.png",
        help="输出图片文件，默认：dy_tip.png",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="只保存图片，不弹出绘图窗口。",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    time, displacement = load_log(input_path)

    minimum_index = int(np.argmin(displacement))

    print(f"读取数据点数：{time.size}")
    print(f"模拟结束时间：{time[-1]:.6g} s")
    print(f"最终端部位移：{displacement[-1]:.6g} m")
    print(
        f"最小端部位移：{displacement[minimum_index]:.6g} m "
        f"（t = {time[minimum_index]:.6g} s）"
    )

    plt.figure(figsize=(9, 6))
    plt.plot(time, displacement, linewidth=1.5, label="Karamelo MPM")

    plt.axhline(0.0, linewidth=0.8)
    plt.scatter(
        time[minimum_index],
        displacement[minimum_index],
        s=35,
        zorder=3,
        label="Minimum displacement",
    )

    plt.xlabel("Time (s)")
    plt.ylabel("Tip displacement (m)")
    plt.title("Cantilever tip displacement")
    plt.xlim(left=0)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"图片已保存：{output_path.resolve()}")

    if not args.quiet:
        plt.show()

    plt.close()


if __name__ == "__main__":
    main()