from enum import IntFlag, StrEnum, auto


class SvdSolverOptions(StrEnum):
    arpack = "arpack"
    lobpcg = "lobpcg"
    auto = "auto"
    randomized = "randomized"

    @classmethod
    def values(cls):
        return [c.value for c in cls]

class ReportingOptions(IntFlag):
    as_progress = auto()
    pdf = auto()

