"""Sequence pairs whose matrix has high scoring cells off the optimal path.

A cell of a Needleman-Wunsch matrix holds the score of the best alignment of
the two *prefixes* that meet there. That is not the same thing as the score
of the best alignment through it: to be on the optimal path a cell also has
to be cheap to *finish* from. A cell can therefore be large and still be a
dead end, and a whole diagonal of such cells can sit next to the traceback
without ever being visited by it.

That is the point these examples are meant to make in a lecture. They are
picked so that the picture drawn by ``alignviz`` shows a clear diagonal
streak of big numbers that the blue traceback arrows walk right past.

Run the module to print them, with the matrix and the streak marked::

    python offpath_examples.py            # the curated examples
    python offpath_examples.py --check    # assert the claims still hold
    python offpath_examples.py --search 20000   # look for new ones

All of them use the defaults of the Streamlit app - global alignment and
``match=3, mismatch=-1, gap=-2`` - so they can be pasted straight into it.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from alignviz import Alignment, Cell, Scoring

__all__ = [
    "EXAMPLES",
    "Example",
    "optimal_cells",
    "best_through",
    "off_path_runs",
    "search",
]

DNA = "ACGT"
DEFAULT_SCORING = Scoring(match=3.0, mismatch=-1.0, gap=-2.0)


# --------------------------------------------------------------------------
# Which cells the optimal alignment actually uses
# --------------------------------------------------------------------------


def optimal_cells(aln: Alignment) -> Set[Cell]:
    """Every cell lying on at least one optimal traceback.

    ``Alignment.tracebacks`` enumerates the paths, of which there can be
    astronomically many; walking the pointer graph backwards from the
    starting cells instead visits each cell once and answers the only
    question asked here, which is whether a cell is used at all.
    """
    seen: Set[Cell] = set()
    stack = list(aln.starts())
    while stack:
        cell = stack.pop()
        if cell in seen:
            continue
        seen.add(cell)
        if aln._stop(cell):
            continue
        for di, dj in aln.trace.get(cell, []):
            stack.append((cell[0] + di, cell[1] + dj))
    return seen


def best_through(aln: Alignment) -> List[List[float]]:
    """``T[i][j]``: the score of the best alignment *through* cell (i, j).

    This is ``S[i][j]`` - the best way in - plus the best way out, so the
    cells with ``T[i][j] == aln.score`` are exactly the optimal ones. The
    gap between a cell's own value and its ``T`` is what makes a tall cell
    a dead end, and ``aln.score - T[i][j]`` says how much an alignment has
    to give up to be routed through it.
    """
    d, m, n = aln.scoring, aln.m, aln.n
    out = [[float("-inf")] * n for _ in range(m)]
    out[m - 1][n - 1] = 0.0
    for i in range(m - 1, -1, -1):
        for j in range(n - 1, -1, -1):
            if i == m - 1 and j == n - 1:
                continue
            best = float("-inf")
            if i + 1 < m and j + 1 < n:
                best = max(best, d(aln.seqA[i], aln.seqB[j]) + out[i + 1][j + 1])
            if i + 1 < m:
                best = max(best, d(aln.seqA[i], "-") + out[i + 1][j])
            if j + 1 < n:
                best = max(best, d("-", aln.seqB[j]) + out[i][j + 1])
            out[i][j] = best
    return [[aln.S[i][j] + out[i][j] for j in range(n)] for i in range(m)]


# --------------------------------------------------------------------------
# Finding the streaks
# --------------------------------------------------------------------------


@dataclass
class Run:
    """A diagonal streak of high cells that the optimal path misses."""

    cells: List[Cell]
    values: List[float]
    separation: int  # how far the streak keeps from the optimal path
    forgone: float  # what routing the alignment through it would cost

    def __len__(self) -> int:
        return len(self.cells)


def off_path_runs(
    aln: Alignment,
    min_len: int = 4,
    above: Optional[float] = None,
    min_separation: int = 2,
) -> List[Run]:
    """Diagonal runs of consecutive off-path cells, longest first.

    Only cells worth at least ``above`` count, which defaults to one point
    more than the optimal score: a streak of cells that each beat the score
    of the best alignment is the version of "high scoring" that needs no
    explaining in a lecture. ``min_separation`` keeps streaks that merely
    graze the traceback out of the way.
    """
    if above is None:
        above = aln.score + 1
    opt = optimal_cells(aln)
    through = best_through(aln)
    runs: List[Run] = []

    def flush(cells: List[Cell]) -> None:
        if len(cells) < min_len:
            return
        sep = min(
            min(abs(i - oi) + abs(j - oj) for oi, oj in opt) for i, j in cells
        )
        if sep < min_separation:
            return
        runs.append(
            Run(
                cells=list(cells),
                values=[aln.S[i][j] for i, j in cells],
                separation=sep,
                forgone=min(aln.score - through[i][j] for i, j in cells),
            )
        )

    for si in range(1, aln.m):
        for sj in range(1, aln.n):
            # start only where the diagonal is broken, so each streak is
            # found once and at its full length
            if si > 1 and sj > 1 and (si - 1, sj - 1) not in opt:
                continue
            i, j, cells = si, sj, []
            while i < aln.m and j < aln.n and (i, j) not in opt:
                if aln.S[i][j] >= above:
                    cells.append((i, j))
                else:
                    flush(cells)
                    cells = []
                i, j = i + 1, j + 1
            flush(cells)

    runs.sort(key=lambda r: (len(r), min(r.values), r.separation), reverse=True)
    return runs


def search(
    trials: int = 20000,
    lengths: Tuple[int, int] = (8, 11),
    scoring: Optional[Scoring] = None,
    min_len: int = 4,
    min_score: float = 3.0,
    seed: Optional[int] = None,
) -> List[Tuple[str, str, Run, Alignment]]:
    """Random DNA pairs that show the effect, best first.

    Sequences of this length keep the drawing readable; ``min_score`` drops
    the pairs whose matrix is all negative numbers, where nothing looks
    high scoring to a student however the arrows run.
    """
    rng = random.Random(seed)
    scoring = scoring or DEFAULT_SCORING
    lo, hi = lengths
    hits: List[Tuple[str, str, Run, Alignment]] = []
    seen: Set[Tuple[str, str]] = set()
    for _ in range(trials):
        a = "".join(rng.choice(DNA) for _ in range(rng.randint(lo, hi)))
        b = "".join(rng.choice(DNA) for _ in range(rng.randint(lo, hi)))
        if (a, b) in seen:
            continue
        seen.add((a, b))
        aln = Alignment(a, b, scoring=scoring, mode="global")
        if aln.score < min_score:
            continue
        runs = off_path_runs(aln, min_len=min_len)
        if runs:
            hits.append((a, b, runs[0], aln))
    hits.sort(
        key=lambda h: (len(h[2]), min(h[2].values) - h[3].score, h[2].separation),
        reverse=True,
    )
    return hits


# --------------------------------------------------------------------------
# The curated examples
# --------------------------------------------------------------------------


@dataclass
class Example:
    seqA: str
    seqB: str
    ridge: List[Cell]
    point: str

    def alignment(self, scoring: Optional[Scoring] = None) -> Alignment:
        return Alignment(
            self.seqA, self.seqB, scoring=scoring or DEFAULT_SCORING, mode="global"
        )


EXAMPLES: List[Example] = [
    Example(
        seqA="CCGCGGTCAC",
        seqB="GCCGTACTTGT",
        ridge=[(6, 3), (7, 4), (8, 5), (9, 6), (10, 7)],
        point=(
            "The biggest number in the matrix is not the answer. The streak "
            "climbs -1, 2, 5, 8, 11 and its 11 is the largest value anywhere, "
            "yet the alignment scores 4, the value in the bottom right corner. "
            "The traceback is a single path and stays four cells clear of the "
            "streak the whole way down."
        ),
    ),
    Example(
        seqA="CTCTTAGTGG",
        seqB="CTTACGCAATC",
        ridge=[(4, 3), (5, 4), (6, 5), (7, 6), (8, 7), (9, 8), (10, 9)],
        point=(
            "Seven cells in a row, 7, 6, 6, 9, 8, 7, 6, every one of them "
            "worth more than the optimal score of 5, and not one of them on "
            "any of the four optimal paths."
        ),
    ),
    Example(
        seqA="CTAGCCCGAT",
        seqB="CGCCCCTCCGC",
        ridge=[(6, 3), (7, 4), (8, 5), (9, 6), (10, 7)],
        point=(
            "The streak 3, 6, 9, 8, 11 runs the length of the lower left half "
            "and peaks at 11 against an optimal score of 5. A good one for "
            "showing that the peak is a strong prefix alignment that is then "
            "far too expensive to finish."
        ),
    ),
    Example(
        seqA="AAGATACAC",
        seqB="CGATAGAAAG",
        ridge=[(2, 1), (3, 2), (4, 3), (5, 4), (6, 5), (7, 6), (8, 7), (9, 8)],
        point=(
            "The streak runs the whole diagonal, -3, 0, 3, 6, 9, 8, 11, 10, "
            "climbing three at a time through the matches while the single "
            "optimal path zigzags around it to a score of 7."
        ),
    ),
]


# --------------------------------------------------------------------------
# Printing
# --------------------------------------------------------------------------


def print_matrix(ex: Example, scoring: Optional[Scoring] = None) -> None:
    """The matrix, with ``*`` on the optimal path and ``#`` on the streak."""
    aln = ex.alignment(scoring)
    opt = optimal_cells(aln)
    ridge = set(ex.ridge)
    a, b = aln.aligned()
    print(f"{ex.seqA} / {ex.seqB}")
    print(f"  score {aln.score:g}, {aln.count_tracebacks()} optimal path(s): {a} / {b}")
    print("      " + "".join(f"{c:>6}" for c in "-" + ex.seqB))
    for i, row in enumerate(aln.S):
        cells = []
        for j, v in enumerate(row):
            mark = "#" if (i, j) in ridge else ("*" if (i, j) in opt else " ")
            cells.append(f"{v:>5g}{mark}")
        print(f"{('-' + ex.seqA)[i]:>6}" + "".join(cells))
    print(f"  # {ex.point}")
    print()


def check(scoring: Optional[Scoring] = None) -> bool:
    """Verify that every curated example still shows what it claims to."""
    ok = True
    for ex in EXAMPLES:
        aln = ex.alignment(scoring)
        opt = optimal_cells(aln)
        stray = [c for c in ex.ridge if c in opt]
        diagonal = all(
            (y[0] - x[0], y[1] - x[1]) == (1, 1)
            for x, y in zip(ex.ridge, ex.ridge[1:])
        )
        peak = max(aln.S[i][j] for i, j in ex.ridge)
        sep = min(
            min(abs(i - oi) + abs(j - oj) for oi, oj in opt) for i, j in ex.ridge
        )
        good = not stray and diagonal and peak > aln.score and sep >= 2
        ok = ok and good
        print(
            f"{'ok  ' if good else 'FAIL'} {ex.seqA} / {ex.seqB}: "
            f"{len(ex.ridge)} cells, peak {peak:g} vs score {aln.score:g}, "
            f"separation {sep}"
            + (f", ON PATH: {stray}" if stray else "")
            + ("" if diagonal else ", NOT DIAGONAL")
        )
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check", action="store_true", help="assert the curated claims still hold"
    )
    parser.add_argument(
        "--search", type=int, metavar="N", help="try N random pairs for new examples"
    )
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    if args.check:
        raise SystemExit(0 if check() else 1)

    if args.search:
        hits = search(trials=args.search, seed=args.seed)
        print(f"{len(hits)} pairs out of {args.search}\n")
        for a, b, run, aln in hits[:15]:
            print(
                f"{a:<12} {b:<12} score {aln.score:>3g}  {len(run)} cells "
                f"{[int(v) for v in run.values]}  separation {run.separation}  "
                f"costs {run.forgone:g} to use"
            )
        return

    for ex in EXAMPLES:
        print_matrix(ex)


if __name__ == "__main__":
    main()
