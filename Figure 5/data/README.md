# Frozen mouse gene-set annotations

- `GO_mouse_annotation.json.gz`: MSigDB M5 GO Biological Process, release 2026.2.Mm.
- `KEGG_mouse_annotation.json.gz`: KEGG mouse pathways, retrieved 2026-11-20.
- `annotation_sources.json`: source URLs, downloaded-content checksums and uncompressed JSON checksums.

These are the exact gene-set snapshots embedded in the revised enrichment notebooks. Each JSON entry maps a term identifier to `[term name, mouse gene symbols]`. The Figure 5 notebook verifies checksums and recomputes enrichment from LANCE8.h5ad without network access. Result tables are written to `../outputs/tables/`.

此处保存修改稿使用的原始冻结注释；运行时校验哈希并重新计算富集结果，无需联网。
