from dsa.binding_interfaces.process_structure.IntervalDict import IntervalDict
from Bio import Align
from dsa.uniprot.uniprot_to_sequence import UniprotToSequence




class SequenceAligner:
    """
    Gap-aware sequence aligner with mismatch/edit fraction threshold.

    The input `sequence` is globally aligned to `reference_sequence`, allowing
    mismatches, insertions, and deletions.

    The best alignment is accepted only if:
      (mismatches + insertions + deletions) <= threshold * len(sequence)

    Mapping is stored as exact-match runs in an IntervalDict:
    for each matched index `i` in sequence, map(i) -> reference_index.
    """

    def __init__(
        self,
        sequence: str,
        uniprot_id: str,
        mismatch_fraction_threshold: float = 0.1,
    ) -> None:
        self.sequence = sequence
        reference_sequence = UniprotToSequence.get_sequence(uniprot_id)
        self.reference_sequence = reference_sequence
        self.mismatch_fraction_threshold = mismatch_fraction_threshold
        self.mapping = IntervalDict()
        self.successful_alignment = False
        self.mismatch_indices = []
        self.mismatch_count: int | None = None
        self.insertion_count: int | None = None
        self.deletion_count: int | None = None

        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = "global"
        self.aligner.match_score = 1.0
        self.aligner.mismatch_score = -1.0
        self.aligner.open_gap_score = -2.0
        self.aligner.extend_gap_score = -0.5
        try:
            self._validate_inputs()
        except Exception as e:
            print(f"Error validating aligner inputs: {e}")
            return None
        self.align()
        # if not self.successful_alignment:
        #     print(f"Alignment failed for {uniprot_id}.")
        #     return None

    def _validate_inputs(self) -> None:
        if not isinstance(self.sequence, str) or not isinstance(self.reference_sequence, str):
            raise TypeError("Both sequence and reference_sequence must be strings.")
        if len(self.sequence) == 0:
            raise ValueError("sequence must not be empty.")
        if len(self.reference_sequence) == 0:
            raise ValueError("reference_sequence must not be empty.")
        if all(aa == 'X' for aa in self.sequence):
            raise ValueError("Not a protein sequence.")
        if all(aa == 'X' for aa in self.reference_sequence):
            raise ValueError("Not a protein sequence.")
        if not (0.0 <= self.mismatch_fraction_threshold <= 1.0):
            raise ValueError("mismatch_fraction_threshold must be between 0.0 and 1.0.")

    @staticmethod
    def sequence_length(sequence: str) -> int:
        return len([aa for aa in sequence if aa != 'X'])

    def _max_allowed_edits(self) -> int:
        return int(self.sequence_length(self.sequence) * self.mismatch_fraction_threshold)

    def _best_alignment(self):
        alignments = self.aligner.align(self.reference_sequence, self.sequence)
        if len(alignments) == 0:
            return None
        return alignments[0]

    def align(self) -> bool:
        """
        Build IntervalDict for exact-match runs in the accepted best alignment.

        Returns:
            True if a valid alignment is found within edit threshold,
            False otherwise.
        """
        alignment = self._best_alignment()
        if alignment is None:
            return False

        self.mapping = IntervalDict()
        self.mismatch_count = 0
        self.insertion_count = 0
        self.deletion_count = 0
        self.edit_count = 0

        coords = alignment.coordinates
        ref_coords = coords[0]
        seq_coords = coords[1]


        for i in range(len(ref_coords) - 1):
            ref_start, ref_end = int(ref_coords[i]), int(ref_coords[i + 1])
            seq_start, seq_end = int(seq_coords[i]), int(seq_coords[i + 1])
            ref_step = ref_end - ref_start
            seq_step = seq_end - seq_start

            if ref_step > 0 and seq_step > 0:
                # Aligned segment (same number of residues on both sides).
                assert ref_step == seq_step
                self.mapping.add(seq_start, seq_end, ref_start-seq_start)
                segment_len = ref_step
                mismatch_indices = [seq_start + i for i in range(segment_len) if self.sequence[seq_start + i] != self.reference_sequence[ref_start + i]]
                self.mismatch_indices.extend(mismatch_indices)
                mismatch_count_without_X = len([i for i in mismatch_indices if self.sequence[i] != 'X'])
                self.mismatch_count += mismatch_count_without_X
            else:
                if ref_step > 0 and seq_step == 0:
                    self.deletion_count += ref_step
                elif seq_step > 0 and ref_step == 0:
                    self.insertion_count += seq_step

    
        if self.mismatch_count > self._max_allowed_edits():
            self.mapping = IntervalDict()
            return False
        else:
            self.successful_alignment = True

        return True

    def map(self, sequence_index: int) -> int | None:
        """
        Return mapped index in reference_sequence for a matching sequence index.
        Returns None for mismatches or if index is outside aligned exact runs.
        """
        mapped_index = self.mapping.map(sequence_index)
        ref_aa = self.reference_sequence[mapped_index]
        seq_aa = self.sequence[sequence_index]
        mismatch = seq_aa != ref_aa
        if mismatch and seq_aa != 'X':
            assert sequence_index in self.mismatch_indices
        return (mapped_index, ref_aa, mismatch)