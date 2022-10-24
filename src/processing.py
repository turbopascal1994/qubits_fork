from argparse import ArgumentParser
from typing import List, Optional


class ResultLog:
    def __init__(
        self,
        cells_number: int,
        sequence: str,
        number_of_cycles: int,
        angle: float,
        fidelity: float,
        execution_time: float,
    ):
        self.cells_number = cells_number
        self.sequence = sequence
        self.number_of_cycles = number_of_cycles
        self.angle = angle
        self.fidelity = fidelity
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
        self.fidelityUpperBound = float(fidelityUpperBound) if fidelityUpperBound is not None else None

    def run(self, array: List[ResultLog]) -> List[ResultLog]:
        ans = []
        for log in array:
            can_take = True
            if (
                self.angle is not None
                and abs(log.angle - self.angle) > self.module
            ):
                can_take = False
            if (
                self.fidelityUpperBound is not None
                and log.fidelity > self.fidelityUpperBound
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
        if len(i) != 6:
            continue
        logs.append(
            ResultLog(int(i[0]), i[1], int(i[2]), float(i[3]), float(i[4]), float(i[5]))
        )
    return logs


def main(args):
    logs = read_file(args.filename)
    result = Filterer(
        angle=args.angle, module=args.module, fidelityUpperBound=args.fidelity
    ).run(logs)
    for log in result:
        print(
            log.cells_number,
            log.number_of_cycles,
            log.angle,
            log.fidelity,
            log.execution_time,
            sep="\t",
        )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--file",
        dest="filename",
        default="result.txt",
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
