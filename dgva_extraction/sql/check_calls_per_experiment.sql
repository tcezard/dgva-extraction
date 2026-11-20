-- dgva-extraction checks, step "calls_per_experiment" (ticket 10, point 5): for one study, the calls with a submitted
-- placement on each assembly, per discovery experiment. With the merged experiment descriptions, it gives the exact
-- number of calls analyses (one per study, submitted assembly and discovery description).
-- Submitted class as in study_records.sql; withdrawn variants left out. Driven by the study index, but it reads every
-- call and placement of the study: run outside busy hours for the large studies (estd214: 6.9 M calls).
-- Binds: study (a study accession).
SELECT
    vs.study_accession,
    p.assembly_accession,
    vce.experiment_id,
    COUNT(DISTINCT vs.variant_accession) AS calls,
    COUNT(*) AS placements
FROM dgva.variant_study vs
    JOIN dgva.variant_call_experiment vce ON vce.variant_call_accession = vs.variant_accession
    JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
WHERE vs.study_accession = :study
    AND vs.region_or_call = 'call'
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
    AND (m.variant_placement_method IS NULL
         OR (m.variant_placement_method NOT LIKE 'Remapped%'
             AND m.variant_placement_method NOT IN ('Remap failed', 'Submitted without placement on assembly',
                                                    'Submitted on non-INSDC assembly')))
GROUP BY vs.study_accession, p.assembly_accession, vce.experiment_id
ORDER BY p.assembly_accession, vce.experiment_id
