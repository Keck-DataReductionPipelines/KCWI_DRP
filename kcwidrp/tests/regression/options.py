"""pytest options of the regression suite (called from the root conftest)."""


def add_options(parser):
    group = parser.getgroup("regression", "KCWI regression suite (T.7)")
    group.addoption(
        "--regression-profile", action="store", default="quick",
        choices=("quick", "full"),
        help="regression entries to run: quick (default) or full (D16)")
    group.addoption(
        "--keep-runs", action="store_true", default=False,
        help="keep redux/ of passing regression runs")
