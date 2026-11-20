-- Ticket 13.3: regions with no placement of their own, whose position is on their supporting calls (translocation
-- regions, nstd83; "Whole variants: measured" in reference/dgva-oracle-schema.md).
-- One row per such region and supporting call placement on the assembly. Only submitted placements of the calls, of
-- non-withdrawn regions and calls. The reader positions the region from these rows.
-- Binds: study (a study accession), assembly (an assembly accession).
SELECT
    vs.variant_accession AS accession,
    rt.sequence_ontology_term AS so_term,
    so.sequence_ontology_id AS so_id,
    r.exchange_variant_region_type AS exchange_type,
    am.assertion_method,
    lv.linked_variant_accession AS call_accession,
    cp.variant_placement_id AS call_placement_id,
    m.variant_placement_method AS placement_method,
    COALESCE(NULLIF(cp.chromosome_name, 'NOT_GIVEN'), ctg.contig_id) AS chromosome,
    COALESCE(st.start_pos, os.outer_start_pos, i_s.inner_start_pos) AS start_pos,
    COALESCE(sp.stop_pos, op.outer_stop_pos, i_p.inner_stop_pos) AS stop_pos
FROM dgva.variant_study vs
    JOIN dgva.variant_region r ON r.variant_region_accession = vs.variant_accession
    LEFT JOIN dgva.variant_region_type_cv rt ON rt.exchange_variant_region_type = r.exchange_variant_region_type
    LEFT JOIN dgva.sequence_ontology_cv so ON so.sequence_ontology_term = rt.sequence_ontology_term
    LEFT JOIN dgva.var_reg_assert_method am ON am.var_reg_assert_method_id = r.var_reg_assert_method_id
    JOIN dgva.linked_variant lv ON lv.variant_region_accession = vs.variant_accession
    JOIN dgva.variant_study cs
        ON cs.variant_accession = lv.linked_variant_accession AND cs.study_accession = vs.study_accession
    JOIN dgva.variant_placement cp ON cp.variant_study_id = cs.variant_study_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_contig ctg ON ctg.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = cp.variant_placement_id
    LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = cp.variant_placement_id
WHERE vs.study_accession = :study
    AND cp.assembly_accession = :assembly
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_placement own WHERE own.variant_study_id = vs.variant_study_id)
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = lv.linked_variant_accession)
    AND rm.variant_placement_id IS NULL
    AND (m.variant_placement_method IS NULL
         OR (m.variant_placement_method LIKE 'Submitted%'
             AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                    'Submitted on non-INSDC assembly')))
ORDER BY vs.variant_accession, chromosome, start_pos, lv.linked_variant_accession, cp.variant_placement_id
