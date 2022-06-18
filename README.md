# evodevo
Mathematical modeling in evolutionary developmental biology.

ark6: Major change.
* The dimension of environmental cues and phenotypes is now equal to nenv (previously nenv + ncells).
* Selection is based on ||p[0:nsel] - env[0:nsel]|| (L1) where nsel < nenv (= ngenes, usually).
* Environmental changes prioritize cue[0:nsel] before cue[nsel:].