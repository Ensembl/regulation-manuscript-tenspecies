# Supporting scripts for Ensembl’s regulatory annotation

This repository contains supporting scripts for the manuscript describing Ensembl's updated regulatory annotation for ten species.

The different scripts and directories are as follows:

| Scripts | Directory | Description |
| --- | --- | ---- |
| Figure 1 | <nobr>[figure_1](figure_1/README.md)</nobr> | Grouping and analysis for Figure&nbsp;1 |
| <nobr>Figures 2&ndash;5</nobr>| <nobr>[figures_2-5](figures_2-5/README.md)</nobr> | Analysis underlying Figures 2&ndash;5 |
| <nobr>Base regulatory annotation</nobr> | <nobr>[regulatory-annotation](regulatory-annotation/README.md)</nobr> | Process primary analysis files to produce promoters, enhancers and unclassified open chromatin regions |
| <nobr>Finalise regulatory annotation</nobr> | <nobr>[ctcf_and_stable_id](ctcf_and_stable_id/README.md)</nobr> | Parse and add CTCF features (where present) to the base annotation, assign stable ids and generate the final GFF3 file. |

#### License

Scripts and code are licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.
