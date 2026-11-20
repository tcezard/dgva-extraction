-- Ticket 10, open point 4 (duplicate IDs): the variants of one study with two or more submitted placements on one
-- assembly, whatever their rank, grouped by type and by what VAR_PLAC_BREAKPOINT_ORDER says about their placements:
-- how many placements have a breakpoint order, and the lowest and highest value; the same for the strand, and how
-- many placements are on the minus strand, or have strand 'not reported'. It tells whether the breakpoint
-- order identifies the variants submitted as two breakpoints (case a), and whether it is absent from the
-- alternative locations (case b) and the pieces (case c).
-- Submitted class as in study_records.sql. Driven by the study index: one study at a time.
-- Binds: study (a study accession).
WITH submitted AS (
    SELECT
        vs.variant_study_id,
        COALESCE(r.exchange_variant_region_type, c.exchange_variant_call_type) AS variant_type,
        p.assembly_accession,
        bo.breakpoint_order,
        str.strand
    FROM dgva.variant_study vs
        JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.var_plac_breakpoint_order bo ON bo.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_strand str ON str.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_region r
            ON vs.region_or_call = 'region' AND r.variant_region_accession = vs.variant_accession
        LEFT JOIN dgva.variant_call c
            ON vs.region_or_call = 'call' AND c.variant_call_accession = vs.variant_accession
    WHERE vs.study_accession = :study
        AND (m.variant_placement_method IS NULL
             OR (m.variant_placement_method NOT LIKE 'Remapped%'
                 AND m.variant_placement_method NOT IN ('Remap failed', 'Submitted without placement on assembly',
                                                        'Submitted on non-INSDC assembly')))
        AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
),
per_variant AS (
    SELECT
        variant_study_id,
        assembly_accession,
        MIN(variant_type) AS variant_type,
        COUNT(*) AS placements,
        COUNT(breakpoint_order) AS with_breakpoint_order,
        MIN(breakpoint_order) AS lowest_breakpoint_order,
        MAX(breakpoint_order) AS highest_breakpoint_order,
        COUNT(strand) AS with_strand,
        SUM(CASE WHEN strand = '-' THEN 1 ELSE 0 END) AS minus_strand,
        SUM(CASE WHEN strand = 'not reported' THEN 1 ELSE 0 END) AS strand_not_reported
    FROM submitted
    GROUP BY variant_study_id, assembly_accession
    HAVING COUNT(*) > 1
)
SELECT
    variant_type,
    placements,
    with_breakpoint_order,
    lowest_breakpoint_order,
    highest_breakpoint_order,
    with_strand,
    minus_strand,
    strand_not_reported,
    COUNT(*) AS variants
FROM per_variant
GROUP BY variant_type, placements, with_breakpoint_order, lowest_breakpoint_order, highest_breakpoint_order,
    with_strand, minus_strand, strand_not_reported
ORDER BY variant_type, placements, with_breakpoint_order, lowest_breakpoint_order, highest_breakpoint_order,
    with_strand, minus_strand, strand_not_reported
