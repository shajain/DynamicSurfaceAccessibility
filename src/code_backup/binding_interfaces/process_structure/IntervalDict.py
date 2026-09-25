from __future__ import annotations
from bisect import bisect_right
from dataclasses import dataclass
import warnings


@dataclass(slots=True)
class _Interval:
    left: int
    right: int
    offset: int   # b_start - a_start


class IntervalDict:
    """
    Maps positions in string A to string B for non-overlapping common substrings.

    Store:  d.add(a_left, a_right, b_left - a_left)
    Query:  d.map(pos_a)  ->  pos_a + offset  (position in B)
    """

    def __init__(self) -> None:
        self._lefts: list[int] = []
        self._intervals: list[_Interval] = []

    def add(self, left: int, right: int, offset: int) -> None:
        if left > right:
            warnings.warn(f"left ({left}) must be <= right ({right})")
            left, right = right, left
        idx = bisect_right(self._lefts, left)
        self._lefts.insert(idx, left)
        self._intervals.insert(idx, _Interval(left, right, offset))

    def __setitem__(self, key: tuple[int, int], offset: int) -> None:
        self.add(*key, offset)

    def map(self, pos_a: int) -> int | None:
        """Return the corresponding position in B, or None if not in any common substring."""
        idx = bisect_right(self._lefts, pos_a) - 1
        if idx < 0:
            return None
        iv = self._intervals[idx]
        if iv.left <= pos_a <= iv.right:
            if iv.offset is None:
                return None
            return pos_a + iv.offset
        #raise ValueError(f"{pos_a} is not in any interval")
        return None

    def __contains__(self, pos_a: int) -> bool:
        return self.map(pos_a) is not None

    def __len__(self) -> int:
        return len(self._intervals)

    def __iter__(self):
        for iv in self._intervals:
            yield iv.left, iv.right, iv.offset

    def __repr__(self) -> str:
        parts = ", ".join(f"[{iv.left},{iv.right}]+{iv.offset}" for iv in self._intervals)
        return f"IntervalDict({{{parts}}})"