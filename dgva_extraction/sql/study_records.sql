-- Ticket 13.2: the records of one study on one assembly.
-- Every placement on the assembly that is not remapped and not in the no-position class, of every non-withdrawn
-- variant of the study. That includes the "conflict" and "other method" placements (none expected): the reader fails
-- the study if one appears, instead of converting it. Placements with no start or no chromosome are returned too; the
-- reader puts them on the side list (W3).
-- Types, SO terms, assertion methods and the assembly are left-joined so that a missing one (W2) reaches the reader,
-- which fails the study, instead of the row silently disappearing as in the legacy exporter.
-- One row per placement and parent region (LINKED_VARIANT): a placement with several parents has several rows. The
-- CIPOS*/CIEND* tables have no primary key, so a duplicated CI row also repeats the placement; the reader checks that
-- repeated rows agree. Rows of one placement are consecutive.
-- Ordered in the database by chromosome name, start, accession and placement, so the writer can stream.
-- Binds: study (a study accession), assembly (an assembly accession).
SELECT
    p.variant_placement_id AS placement_id,
    vs.variant_accession AS accession,
    vs.region_or_call,
    COALESCE(rt.sequence_ontology_term, ct_type.sequence_ontology_term) AS so_term,
    so.sequence_ontology_id AS so_id,
    COALESCE(r.exchange_variant_region_type, c.exchange_variant_call_type) AS exchange_type,
    CASE WHEN r.variant_region_accession IS NULL AND c.variant_call_accession IS NULL THEN 0 ELSE 1 END
        AS in_type_table,
    am.assertion_method,
    CASE WHEN ga.assembly_accession IS NULL THEN 0 ELSE 1 END AS assembly_known,
    m.variant_placement_method AS placement_method,
    CASE WHEN rm.variant_placement_id IS NULL THEN 0 ELSE 1 END AS has_remap_score,
    NULLIF(p.chromosome_name, 'NOT_GIVEN') AS chromosome_name,
    ctg.contig_id,
    COALESCE(NULLIF(p.chromosome_name, 'NOT_GIVEN'), ctg.contig_id) AS chromosome,
    st.start_pos,
    sp.stop_pos,
    os.outer_start_pos,
    i_s.inner_start_pos,
    i_p.inner_stop_pos,
    op.outer_stop_pos,
    cs.ci_minus AS ci_start_minus,
    cs.ci_plus AS ci_start_plus,
    ce.ci_minus AS ci_stop_minus,
    ce.ci_plus AS ci_stop_plus,
    pl.ciposleft,
    pr.ciposright,
    el.ciendleft,
    er.ciendright,
    lv.variant_region_accession AS parent_accession,
    CASE WHEN ps.variant_study_id IS NULL THEN 0 ELSE 1 END AS parent_in_study,
    CASE WHEN pw.variant_accession IS NULL THEN 0 ELSE 1 END AS parent_withdrawn
FROM dgva.variant_study vs
    JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = p.assembly_accession
    LEFT JOIN dgva.variant_region r
        ON vs.region_or_call = 'region' AND r.variant_region_accession = vs.variant_accession
    LEFT JOIN dgva.variant_call c
        ON vs.region_or_call = 'call' AND c.variant_call_accession = vs.variant_accession
    LEFT JOIN dgva.variant_region_type_cv rt ON rt.exchange_variant_region_type = r.exchange_variant_region_type
    LEFT JOIN dgva.variant_call_type_cv ct_type ON ct_type.exchange_variant_call_type = c.exchange_variant_call_type
    LEFT JOIN dgva.sequence_ontology_cv so
        ON so.sequence_ontology_term = COALESCE(rt.sequence_ontology_term, ct_type.sequence_ontology_term)
    LEFT JOIN dgva.var_reg_assert_method am ON am.var_reg_assert_method_id = r.var_reg_assert_method_id
    LEFT JOIN dgva.variant_placement_contig ctg ON ctg.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_start cs ON cs.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ci_stop ce ON ce.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciposleft pl ON pl.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciposright pr ON pr.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciendleft el ON el.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_ciendright er ON er.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.linked_variant lv ON lv.linked_variant_accession = vs.variant_accession
    LEFT JOIN dgva.variant_study ps
        ON ps.variant_accession = lv.variant_region_accession AND ps.study_accession = vs.study_accession
    LEFT JOIN dgva.variant_withdrawn pw ON pw.variant_accession = lv.variant_region_accession
WHERE vs.study_accession = :study
    AND p.assembly_accession = :assembly
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    AND (m.variant_placement_method IS NULL
         OR (m.variant_placement_method NOT LIKE 'Remapped%'
             AND m.variant_placement_method NOT IN ('Remap failed', 'Submitted without placement on assembly',
                                                    'Submitted on non-INSDC assembly')))
ORDER BY
    COALESCE(NULLIF(p.chromosome_name, 'NOT_GIVEN'), ctg.contig_id),
    COALESCE(st.start_pos, os.outer_start_pos, i_s.inner_start_pos),
    vs.variant_accession,
    p.variant_placement_id,
    lv.variant_region_accession
