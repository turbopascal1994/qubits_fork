from argparse import ArgumentParser
from typing import List, Optional


class ResultLog:
    def __init__(
        self,
        cells_number: int,
        sequence: str,
        number_of_cycles: int,
        angle: float,
        fidelity1: float,
        fidelity2: float,
        fidelity3: float,
        fidelity12: float,
        execution_time: float,
    ):
        self.cells_number = cells_number
        self.sequence = sequence
        self.number_of_cycles = number_of_cycles
        self.angle = angle
        self.fidelity1 = fidelity1
        self.fidelity2 = fidelity2
        self.fidelity3 = fidelity3
        self.fidelity12 = fidelity12
        self.execution_time = execution_time


class Filterer:
    def __init__(
        self,
        angle: Optional[str] = None,
        module: Optional[str] = None,
        fidelityUpperBound: Optional[str] = None,
    ):
        self.angle = float(angle) if angle is not None else None
        self.module = float(module) if module is not None else None
        self.fidelityUpperBound = (
            float(fidelityUpperBound) if fidelityUpperBound is not None else None
        )

    def run(self, array: List[ResultLog]) -> List[ResultLog]:
        ans = []
        for log in array:
            can_take = True
            if self.angle is not None and abs(log.angle - self.angle) > self.module:
                can_take = False
            if (
                self.fidelityUpperBound is not None
                and log.fidelity3 > self.fidelityUpperBound
            ):
                can_take = False
            if can_take:
                ans.append(log)
        return ans


def read_file(filename):
    with open(filename, "r") as f:
        input = f.read()

    input = input.split("\n")
    for index, i in enumerate(input):
        input[index] = i.split()
    logs = []
    for i in input:
        if len(i) != 9:
            continue
        logs.append(
            ResultLog(
                cells_number=int(i[0]),
                sequence=i[1],
                number_of_cycles=int(i[2]),
                angle=float(i[3]),
                fidelity1=float(i[4]),
                fidelity2=float(i[5]),
                fidelity3=float(i[6]),
                fidelity12=float(i[7]),
                execution_time=float(i[8]),
            )
        )
    return logs


class ResultLogFilename:
    def __init__(self, filename):
        filename = filename[:-4]
        args = filename.split("_")
        self.params = {}
        for i in args:
            self.params[i.split("=")[0]] = i.split("=")[1]


def main(args):
    import os
    from collections import defaultdict

    files = os.listdir(args.folder)
    final_dict = defaultdict(list)
    for file in files:
        filename_params = ResultLogFilename(file)
        logs = read_file(os.path.join(args.folder, file))
        filtered = Filterer(
            angle=args.angle, module=args.module, fidelityUpperBound=args.fidelity
        ).run(logs)
        filtered = sorted(filtered, key=lambda x: x.fidelity3)
        if len(filtered) > 0:
            final_dict[filename_params.params["w01"]].append(filtered[0])

    for key, val in final_dict.items():
        with open(key + ".txt", "w") as f:
            for log in val:
                f.write(
                    str(log.cells_number)
                    + "\t"
                    + str(log.sequence)
                    + "\t"
                    + str(log.number_of_cycles)
                    + "\t"
                    + str(log.angle)
                    + "\t"
                    + str(log.fidelity1)
                    + "\t"
                    + str(log.fidelity2)
                    + "\t"
                    + str(log.fidelity3)
                    + "\t"
                    + str(log.fidelity12)
                    + "\t"
                    + str(log.execution_time)
                    + "\n"
                )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--folder",
        dest="folder",
        help="Path to the file to be processed.",
    )
    parser.add_argument(
        "--angle",
        dest="angle",
        default=None,
        help="Required angle. If None, then filtering by this attribute is not necessary.",
    )
    parser.add_argument(
        "--module",
        dest="module",
        default=None,
        help="The module that the filtering will be relative to. After filtering, the values will remain such that abs(x - angle) <= module. If angle is None, then skip this attribute.",
    )
    parser.add_argument(
        "--fidelity",
        dest="fidelity",
        default=None,
        help="Upper limit of fidelity. After filtering, fidelity values that are less or equal than the value will remain. If None, then filtering by this attribute is not necessary.",
    )
    args = parser.parse_args()

    main(args)
