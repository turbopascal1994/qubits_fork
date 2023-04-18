import cmath
from argparse import ArgumentParser

import qutip
from tqdm import tqdm


def transform_to_coordinates(wf0: complex, wf1: complex):
    polar0, polar1 = cmath.polar(wf0), cmath.polar(wf1)
    phi = polar1[1] - polar0[1]
    if phi < 0:
        phi += 2 * cmath.pi
    assert 0 <= phi < 2 * cmath.pi, phi
    teta = 2 * cmath.acos(polar1[0]).real
    assert 0 <= teta <= cmath.pi, teta

    return (
        cmath.sin(teta) * cmath.cos(phi),
        cmath.sin(teta) * cmath.sin(phi),
        cmath.cos(teta),
    )


def read_file(filename):
    with open(filename, "r") as f:
        input = f.read()

    input = input.split("\n")
    x = [[] for _ in range(6)]
    y = [[] for _ in range(6)]
    z = [[] for _ in range(6)]
    color = []
    for index, i in enumerate(input):
        input[index] = list(map(float, i.split()))
        color.append(int(input[index][-1]))
        X, Y, Z = transform_to_coordinates(
            complex(input[index][0], input[index][1]),
            complex(input[index][2], input[index][3]),
        )
        x[color[-1]].append(X)
        y[color[-1]].append(Y)
        z[color[-1]].append(Z)
    xx = []
    yy = []
    zz = []
    for i in range(len(z[-1])):
        for j in range(6):
            xx.append(x[j][i])
            yy.append(y[j][i])
            zz.append(z[j][i])
    return xx, yy, zz


def get_alpha(start, end, cur_step, steps):
    return start + (end - start) / steps * cur_step


class Buffer:
    def __init__(self, size: int):
        self.size = size
        self._buffer = []

    def add_points(self, x: list, y: list, z: list):
        self._buffer.append([x, y, z])
        if len(self._buffer) > self.size:
            self._buffer = self._buffer[1:]

    def visualize(self, sphere: qutip.Bloch):
        for index, i in enumerate(self._buffer):
            sphere.add_points(
                i, meth="m", alpha=get_alpha(0.3, 1.0, index, len(self._buffer))
            )


def main(args):
    x, y, z = read_file(args.path)
    colors = ["red", "green", "blue", "cyan", "magenta", "yellow"]
    sphere = qutip.Bloch()
    # sphere.view = [-40, 30]
    sphere.point_color = colors
    sphere.vector_color = colors
    sphere.point_marker = ["o"]
    buffer = Buffer(size=10)
    for i in tqdm(range(0, len(x) // 6)):
        sphere.clear()
        buffer.add_points(
            x[6 * i : 6 * (i + 1)],
            y[6 * i : 6 * (i + 1)],
            z[6 * i : 6 * (i + 1)],
        )
        buffer.visualize(sphere)
        sphere.save(
            dirc="temp"
        )  # saving images to temp directory in current working directory


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--path",
        dest="path",
        help="Path to the file to be processed.",
    )
    args = parser.parse_args()
    main(args)
"""
python visualize.py --path D:\Diplom\ohfss\src\tmp.txt
ffmpeg -i temp/bloch_%01d.png -pix_fmt yuv420p -r 40 bloch.mp4
"""
