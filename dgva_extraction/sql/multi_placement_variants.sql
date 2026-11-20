-- Ticket 10, point 4: the variants of one study with two or more submitted placements on one assembly, with what
-- decides how they are written: their type, coordinates, placement rank, breakpoint order (From/To) and strand.
-- The reader reads this before study_records.sql (a few thousand rows at most: 58,284 variants in nstd115) and plans
-- each variant: a breakend pair, single breakends, one record spanning two breakpoints, or records linked by EVENT.
-- Same placements as study_records.sql: submitted, of non-withdrawn variants, with a start and a chromosome (W3
-- placements are on the side list and don't count). Start and stop are taken as the reader takes them: start, then
-- outer, then inner.
-- Driven by the study index, like study_records.sql.
-- Binds: study (a study accession), assembly (an assembly accession).
WITH placed AS (
    SELECT
        vs.variant_accession AS accession,
        p.variant_placement_id AS placement_id,
        COALESCE(rt.sequence_ontology_term, ct_type.sequence_ontology_term) AS so_term,
        COALESCE(NULLIF(p.chromosome_name, 'NOT_GIVEN'), ctg.contig_id) AS chromosome,
        COALESCE(st.start_pos, os.outer_start_pos, i_s.inner_start_pos) AS start_pos,
        COALESCE(sp.stop_pos, op.outer_stop_pos, i_p.inner_stop_pos) AS stop_pos,
        rk.placement_rank,
        bo.breakpoint_order,
        str.strand
    FROM dgva.variant_study vs
        JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_region r
            ON vs.region_or_call = 'region' AND r.variant_region_accession = vs.variant_accession
        LEFT JOIN dgva.variant_call c
            ON vs.region_or_call = 'call' AND c.variant_call_accession = vs.variant_accession
        LEFT JOIN dgva.variant_region_type_cv rt ON rt.exchange_variant_region_type = r.exchange_variant_region_type
        LEFT JOIN dgva.variant_call_type_cv ct_type ON ct_type.exchange_variant_call_type = c.exchange_variant_call_type
        LEFT JOIN dgva.variant_placement_contig ctg ON ctg.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_stop sp ON sp.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_inner_stop i_p ON i_p.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_outer_stop op ON op.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.var_plac_placement_rank rk ON rk.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.var_plac_breakpoint_order bo ON bo.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_strand str ON str.variant_placement_id = p.variant_placement_id
    WHERE vs.study_accession = :study
        AND p.assembly_accession = :assembly
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method NOT LIKE 'Remapped%'
                 AND m.variant_placement_method NOT IN ('Remap failed', 'Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
),
positioned AS (
    SELECT * FROM placed WHERE chromosome IS NOT NULL AND start_pos IS NOT NULL
),
several AS (
    SELECT accession FROM positioned GROUP BY accession HAVING COUNT(*) > 1
)
SELECT positioned.*
FROM positioned
    JOIN several ON several.accession = positioned.accession
ORDER BY positioned.accession, positioned.placement_id
