Utility scripts
================
- ``blast_to_raw.py``: a powerful BLAST filter, can filter tandem matches,
  repetitive matches, and weak matches. often run before ``quota_alignment.py``,
  however can also run on its own (and writes a filtered BLAST result).
- ``qa_plot.py``: after running ``quota_alignment.py``, visualize the dot plot
- ``qa_to_pairs.py``: after running ``quota_alignment.py``, generate the anchor
  pairs
- ``synteny_score.py``: call "syntelog" versus "gray" genes
  (transposed genes).
- ``synteny_liftover.py``: extract hits in the blast_file that are in the vicinity of the anchors. Useful to improve your anchor list (so that it is more comprehensive). Typical use include an anchor list defined by all-vs-all protein BLAST, then add the noncoding stuff (like miRNAs) to the list.
