-- Ticket 10, open point 4 (duplicate IDs): every placement of one variant of one study, submitted or remapped, with
-- the variant's type, and the fields that tell the placements apart (rank, placements per assembly, assembly unit,
-- alt status, breakpoint order, strand).
-- Binds: variant_study_id (VARIANT_STUDY.VARIANT_STUDY_ID).
SELECT
    p.variant_placement_id AS placement_id,
    COALESCE(r.exchange_variant_region_type, c.exchange_variant_call_type) AS variant_type,
    p.assembly_accession,
    m.variant_placement_method AS placement_method,
    CASE WHEN rm.variant_placement_id IS NULL THEN 0 ELSE 1 END AS has_remap_score,
    rk.placement_rank,
    pa.placements_per_assembly,
    au.assembly_unit,
    alt.alt_status,
    bo.breakpoint_order,
    str.strand,
    COALESCE(NULLIF(p.chromosome_name, 'NOT_GIVEN'), ctg.contig_id) AS chromosome,
    st.start_pos,
    sp.stop_pos,
    os.outer_start_pos,
    i_s.inner_start_pos,
    i_p.inner_stop_pos,
    op.outer_stop_pos
FROM dgva.variant_placement p
    JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
    LEFT JOIN dgva.variant_region r
        ON vs.region_or_call = 'region' AND r.variant_region_accession = vs.variant_accession
    LEFT JOIN dgva.variant_call c
        ON vs.region_or_call = 'call' AND c.variant_call_accession = vs.variant_accession
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.var_plac_placement_rank rk ON rk.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.var_plac_placs_per_assem pa ON pa.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_assem_unit au ON au.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_alt_status alt ON alt.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.var_plac_breakpoint_order bo ON bo.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_strand str ON str.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_contig ctg ON ctg.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = p.variant_placement_id
WHERE p.variant_study_id = :variant_study_id
ORDER BY p.assembly_accession, p.variant_placement_id
