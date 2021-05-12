

# Store WT sequence...
# How fast can I read in fastq files?
#

def count_variant(self, variant_dna, include_indels=True):
    """
    Identifies mutations and counts the *variant_dna* sequence.
    The algorithm attempts to call variants by comparing base-by-base.
    If the *variant_dna* and wild type DNA are different lengths, or if
    there are an excess of mismatches (indicating a possible indel), local
    alignment is performed using :py:meth:`align_variant` if this option
    has been selected in the configuration.

    Each variant is stored as a tab-delimited string of mutations in HGVS_
    format. Returns a list of HGVS_ variant strings. Returns an empty list
    if the variant is wild type. Returns None if the variant was discarded
    due to excess mismatches.
    """
    if not re.match("^[ACGTNXacgtnx]+$", variant_dna):
        raise ValueError(
            "Variant DNA sequence contains unexpected "
            "characters [{}]".format(self.name)
        )

    variant_dna = variant_dna.upper()

    if len(variant_dna) != len(self.wt.dna_seq):
        if self.aligner is not None:
            mutations = self.align_variant(variant_dna)
        else:
            return None
    else:
        mutations = list()
        for i in xrange(len(variant_dna)):
            if variant_dna[i] != self.wt.dna_seq[i]:
                mutations.append(
                    (
                        i,
                        "{pre}>{post}".format(
                            pre=self.wt.dna_seq[i], post=variant_dna[i]
                        ),
                    )
                )
                if len(mutations) > self.max_mutations:
                    if self.aligner is not None:
                        mutations = self.align_variant(variant_dna)
                        if len(mutations) > self.max_mutations:
                            # too many mutations post-alignment
                            return None
                        else:
                            # stop looping over this variant
                            break
                    else:
                        # too many mutations and not using aligner
                        return None

    mutation_strings = list()
    if self.is_coding():
        variant_protein = ""
        for i in xrange(0, len(variant_dna), 3):
            try:
                variant_protein += CODON_TABLE[variant_dna[i : i + 3]]
            except KeyError:  # garbage codon due to indel, X, or N
                variant_protein += "?"

        for pos, change in mutations:
            ref_dna_pos = pos + self.wt.dna_offset + 1
            ref_pro_pos = pos / 3 + self.wt.protein_offset + 1
            mut = "c.{pos}{change}".format(pos=ref_dna_pos, change=change)
            if has_indel(change):
                mut += " (p.{pre}{pos}fs)".format(
                    pre=AA_CODES[self.wt.protein_seq[pos / 3]], pos=ref_pro_pos
                )
            elif variant_protein[pos / 3] == self.wt.protein_seq[pos / 3]:
                mut += " (p.=)"
            else:
                mut += " (p.{pre}{pos}{post})".format(
                    pre=AA_CODES[self.wt.protein_seq[pos / 3]],
                    pos=ref_pro_pos,
                    post=AA_CODES[variant_protein[pos / 3]],
                )
            mutation_strings.append(mut)
    else:
        for pos, change in mutations:
            ref_dna_pos = pos + self.wt.dna_offset + 1
            mut = "n.{pos}{change}".format(pos=ref_dna_pos, change=change)
            mutation_strings.append(mut)

    if len(mutation_strings) > 0:
        variant_string = ", ".join(mutation_strings)
    else:
        variant_string = WILD_TYPE_VARIANT
    return variant_string

def count_synonymous(self):
    """
    Combine counts for synonymous variants (defined as variants that differ
    at the nucleotide level but not at the amino acid level) and store them
    under the ``synonymous`` label.

    This method should be called only after ``variants`` have been counted.

    .. note:: The total number of ``synonymous`` variants may be greater \
    than the total number of ``variants`` after filtering. This is \
    because ``variants`` are combined into ``synonymous`` entries at the \
    :py:class:`~seqlib.seqlib.SeqLib` level before count-based filtering, \
    allowing filtered ``variants`` to contribute counts to their \
    ``synonymous`` entry.
    """
    if not self.is_coding():
        self.logger.warning("Cannot count synonymous mutations in noncoding data")
        return

    if self.check_store("/main/synonymous/counts"):
        return

    self.logger.info("Counting synonymous variants")

    df_dict = dict()

    for variant, count in self.store["/main/variants/counts"].iterrows():
        if variant == WILD_TYPE_VARIANT:
            df_dict[variant] = count["count"]
        else:
            variant = protein_variant(variant)
            if len(variant) == 0:
                variant = SYNONYMOUS_VARIANT
            try:
                df_dict[variant] += count["count"]
            except KeyError:
                df_dict[variant] = count["count"]

    self.save_counts("synonymous", df_dict, raw=False)
    del df_dict
