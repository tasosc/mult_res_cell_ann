from enum import StrEnum


class SvdSolverOptions(StrEnum):
    arpack = "arpack"
    lobpcg = "lobpcg"
    auto = "auto"
    randomized = "randomized"

    @classmethod
    def values(cls):
        return [c.value for c in cls]
