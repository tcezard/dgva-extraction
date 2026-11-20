-- Ticket 13.3: the variants of one study that are in no VCF, for the side list (sprint 4 turns it into the DGVa JSON).
--   withdrawn:          the accession is in VARIANT_WITHDRAWN (with its reason and the study that withdrew it)
--   remapped_only:      no submitted placement, only DGVa remappings (nstd35)
--   no_position_only:   only 'Remap failed' or 'Submitted without placement ...' placements
--   no_placement:       no placement at all (estd180; nstd83 translocation regions, which the reader may still
--                       position from their calls)
-- Placement-level entries (a submitted placement with no start or no chromosome, duplicate placements) are added by
-- the reader. Withdrawn variants are listed whatever placements they have.
-- Binds: study (a study accession).
SELECT
    accession,
    region_or_call,
    CASE
        WHEN withdrawal_reason IS NOT NULL THEN 'withdrawn'
        WHEN n_placements = 0 THEN 'no_placement'
        WHEN n_remapped > 0 THEN 'remapped_only'
        ELSE 'no_position_only'
    END AS reason,
    withdrawal_reason,
    withdrawn_by_study,
    n_placements,
    n_remapped,
    n_no_position
FROM (
    SELECT
        vs.variant_accession AS accession,
        vs.region_or_call,
        w.withdrawal_reason,
        w.study_accession AS withdrawn_by_study,
        COUNT(p.variant_placement_id) AS n_placements,
        SUM(CASE WHEN p.variant_placement_id IS NOT NULL AND rm.variant_placement_id IS NULL
                      AND (m.variant_placement_method IS NULL
                           OR (m.variant_placement_method LIKE 'Submitted%'
                               AND m.variant_placement_method NOT IN ('Submitted without placement on assembly',
                                                                      'Submitted on non-INSDC assembly')))
                 THEN 1 ELSE 0 END) AS n_submitted,
        SUM(CASE WHEN m.variant_placement_method LIKE 'Remapped%' THEN 1 ELSE 0 END) AS n_remapped,
        SUM(CASE WHEN m.variant_placement_method IN ('Remap failed', 'Submitted without placement on assembly',
                                                     'Submitted on non-INSDC assembly') THEN 1 ELSE 0 END)
            AS n_no_position
    FROM dgva.variant_study vs
        LEFT JOIN dgva.variant_withdrawn w ON w.variant_accession = vs.variant_accession
        LEFT JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
        LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
        LEFT JOIN dgva.variant_placement_remap rm ON rm.variant_placement_id = p.variant_placement_id
    WHERE vs.study_accession = :study
    GROUP BY vs.variant_accession, vs.region_or_call, w.withdrawal_reason, w.study_accession
) variant_placements
WHERE withdrawal_reason IS NOT NULL OR n_submitted = 0
ORDER BY accession
