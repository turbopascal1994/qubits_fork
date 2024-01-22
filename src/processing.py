from argparse import ArgumentParser
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class ResultLog:
    cells_number: int
    sequence: str
    number_of_cycles: int
    iterations: int
    fidelity: float
    leak: float
    execution_time: float


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
            # if self.angle is not None and abs(log.angle - self.angle) > self.module:
            #     can_take = False
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
        if len(i) != 7:
            continue
        logs.append(
            ResultLog(
                cells_number=int(i[0]),
                sequence=i[1],
                number_of_cycles=int(i[2]),
                iterations=int(i[3]),
                fidelity=float(i[4]),
                leak=float(i[5]),
                execution_time=float(i[6]),
            )
        )
    return logs


class ResultLogFilename:
    def __init__(self, filename):
        filename = filename[:-4]
        args = filename.split("_")
        self.params = {}
        for i in args:
            splited = i.split('=')
            self.params[splited[0]] = splited[1]


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
        filtered = sorted(filtered, key=lambda x: x.fidelity)
        if len(filtered) > 0:
            filename = f"N={int(filename_params.params['RL'])}, M={int(filename_params.params['L'])}, w01={filename_params.params['w01']}, w12={filename_params.params['w12']}, wt={filename_params.params['wt']}, angle={filename_params.params['Angle']}"
            # filename = f"{filename_params.params['w01']}"
            final_dict[filename].append(filtered[0])
            for i in range(1, len(filtered)):
                # if filtered[i].number_of_cycles == filtered[i - 1].number_of_cycles:
                #     continue
                final_dict[filename].append(filtered[i])

    for key, val in final_dict.items():
        val = sorted(val, key=lambda x: x.cells_number)
        with open(key + ".txt", "w") as f:
            for log in val:
                log: ResultLog
                f.write(
                    str(log.cells_number)
                    + "\t"
                    + str(log.sequence)
                    + "\t"
                    + str(log.number_of_cycles)
                    + "\t"
                    + str(log.iterations)
                    + "\t"
                    + str(log.fidelity)
                    + "\t"
                    + str(log.leak)
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
        default=0.0001,
        help="Upper limit of fidelity. After filtering, leak values that are less or equal than the value will remain. If None, then filtering by this attribute is not necessary.",
    )
    args = parser.parse_args()

    main(args)
