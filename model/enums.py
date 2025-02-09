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

class Activity(StrEnum):
    NONE = "",
    PARSE_DATASET = "parse_dataset",
    UPLOAD_DATASET = "upload_dataset",
    PP_QC = "pp_qc",
    PP_NORM = "pp_nrom",
    PP_FEATURE = "pp_feature",
    PP_REDUCTION = "pp_reduction",
    PP_VISUALIAZTION = "pp_visualiaztion",
    SI_CLUSTERING = "si_clustering",
    SI_ANNOTATION = "si_annotation",
    SI_FILE = "si_file",
    END = auto()
