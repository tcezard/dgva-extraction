-- Ticket 13.1: the study and assembly list, the list of what to convert (plan step 3.3).
-- One row per study and assembly holding placements of that study; a study with no placement (estd235) gets one row
-- with a null assembly. Placements are classed by the ticket 2 rule (plan step 1.4):
--   no_position: method 'Remap failed', 'Submitted without placement on assembly', 'Submitted on non-INSDC assembly'
--   remapped:    method starts with 'Remapped'
--   conflict:    method null or 'Submitted ...', with a remap score (none expected; reported, never converted)
--   submitted:   method null or 'Submitted ...'
--   other:       any other method (none expected; reported, never converted)
-- Withdrawn variants (any accession in VARIANT_WITHDRAWN, as the legacy exporters do) are counted apart and left out of
-- every other count. An assembly is a submitted assembly, to convert, when submitted_placements > 0.
-- Check: the sum of submitted_variants over all rows is 44,156,767 (sprint 2, ticket 6); no variant is submitted on
-- two assemblies.
-- Binds: study (a study accession, or ALL). With ALL it reads all 96.6 M placements: run outside busy hours.
WITH placement_class AS (
    SELECT
        p.variant_placement_id,
        p.variant_study_id,
        p.assembly_accession,
        m.variant_placement_method AS placement_method,
        CASE
            WHEN m.variant_placement_method IN ('Remap failed', 'Submitted without placement on assembly',
                                                'Submitted on non-INSDC assembly') THEN 'no_position'
            WHEN m.variant_placement_method LIKE 'Remapped%' THEN 'remapped'
            WHEN (m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%')
                 AND r.variant_placement_id IS NOT NULL THEN 'conflict'
            WHEN m.variant_placement_method IS NULL OR m.variant_placement_method LIKE 'Submitted%' THEN 'submitted'
            ELSE 'other'
        END AS placement_class,
        CASE WHEN COALESCE(st.start_pos, os.outer_start_pos, i_s.inner_start_pos) IS NULL THEN 0 ELSE 1 END
            AS has_start,
        CASE WHEN NULLIF(p.chromosome_name, 'NOT_GIVEN') IS NOT NULL OR ct.contig_id IS NOT NULL THEN 1 ELSE 0 END
            AS has_name
    FROM dgva.variant_placement p
        JOIN dgva.variant_study vs ON vs.variant_study_id = p.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap r ON r.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_start st ON st.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_outer_start os ON os.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_inner_start i_s ON i_s.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_contig ct ON ct.variant_placement_id = p.variant_placement_id
    WHERE :study = 'ALL' OR vs.study_accession = :study
),
variant AS (
    SELECT
        vs.variant_study_id,
        vs.study_accession,
        CASE WHEN w.variant_accession IS NULL THEN 0 ELSE 1 END AS withdrawn
    FROM dgva.variant_study vs
        LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = vs.variant_accession
    WHERE :study = 'ALL' OR vs.study_accession = :study
),
study_assembly AS (
    SELECT
        v.study_accession,
        pc.assembly_accession,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'submitted' THEN 1 ELSE 0 END)
            AS submitted_placements,
        COUNT(DISTINCT CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'submitted'
                            THEN v.variant_study_id END) AS submitted_variants,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'submitted' AND pc.placement_method IS NULL
                 THEN 1 ELSE 0 END) AS submitted_null_method,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'submitted' AND pc.has_start = 0
                 THEN 1 ELSE 0 END) AS submitted_without_start,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'submitted' AND pc.has_start = 1
                      AND pc.has_name = 0 THEN 1 ELSE 0 END) AS submitted_without_chromosome,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'remapped' THEN 1 ELSE 0 END) AS remapped_placements,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'no_position' THEN 1 ELSE 0 END)
            AS no_position_placements,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'conflict' THEN 1 ELSE 0 END) AS conflict_placements,
        SUM(CASE WHEN v.withdrawn = 0 AND pc.placement_class = 'other' THEN 1 ELSE 0 END) AS other_method_placements,
        SUM(v.withdrawn) AS withdrawn_placements
    FROM variant v
        JOIN placement_class pc ON pc.variant_study_id = v.variant_study_id
    GROUP BY v.study_accession, pc.assembly_accession
),
variant_summary AS (
    SELECT
        v.study_accession,
        v.variant_study_id,
        v.withdrawn,
        COUNT(pc.variant_placement_id) AS n_placements,
        SUM(CASE WHEN pc.placement_class = 'submitted' THEN 1 ELSE 0 END) AS n_submitted,
        SUM(CASE WHEN pc.placement_class = 'remapped' THEN 1 ELSE 0 END) AS n_remapped
    FROM variant v
        LEFT JOIN placement_class pc ON pc.variant_study_id = v.variant_study_id
    GROUP BY v.study_accession, v.variant_study_id, v.withdrawn
),
study_summary AS (
    SELECT
        study_accession,
        COUNT(*) AS study_variants,
        SUM(withdrawn) AS study_withdrawn_variants,
        SUM(CASE WHEN withdrawn = 0 AND n_submitted > 0 THEN 1 ELSE 0 END) AS study_submitted_variants,
        SUM(CASE WHEN withdrawn = 0 AND n_submitted = 0 AND n_remapped > 0 THEN 1 ELSE 0 END)
            AS study_remapped_only_variants,
        SUM(CASE WHEN withdrawn = 0 AND n_placements > 0 AND n_submitted = 0 AND n_remapped = 0 THEN 1 ELSE 0 END)
            AS study_no_position_variants,
        SUM(CASE WHEN withdrawn = 0 AND n_placements = 0 THEN 1 ELSE 0 END) AS study_unplaced_variants
    FROM variant_summary
    GROUP BY study_accession
)
SELECT
    s.study_accession,
    s.display_name,
    s.study_type,
    s.hold_date,
    s.study_complete,
    ss.study_variants,
    ss.study_withdrawn_variants,
    ss.study_submitted_variants,
    ss.study_remapped_only_variants,
    ss.study_no_position_variants,
    ss.study_unplaced_variants,
    sa.assembly_accession,
    ga.assembly_name,
    ga.assembly_taxonomy_id,
    CASE WHEN sa.assembly_accession IS NOT NULL AND ga.assembly_accession IS NULL THEN 1 ELSE 0 END
        AS assembly_unknown,
    sa.submitted_placements,
    sa.submitted_variants,
    sa.submitted_null_method,
    sa.submitted_without_start,
    sa.submitted_without_chromosome,
    sa.remapped_placements,
    sa.no_position_placements,
    sa.conflict_placements,
    sa.other_method_placements,
    sa.withdrawn_placements
FROM dgva.dgva_study s
    LEFT JOIN study_summary ss ON ss.study_accession = s.study_accession
    LEFT JOIN study_assembly sa ON sa.study_accession = s.study_accession
    LEFT JOIN dgva.genome_assembly ga ON ga.assembly_accession = sa.assembly_accession
WHERE :study = 'ALL' OR s.study_accession = :study
ORDER BY s.study_accession, sa.assembly_accession
