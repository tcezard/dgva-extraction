-- Ticket 10, open point 4 (duplicate IDs): variants with a submitted placement ranked NotBest, i.e. a variant that
-- DGVa also placed elsewhere on an assembly it was submitted on. dgva-extraction duplicate-placements reads a few of
-- these and stops; it then reads every placement of each one with variant_placements.sql.
-- One study at a time, so that the study index drives the query: over all studies it reads the whole rank table
-- (about 3 minutes; 10 studies have such variants: see docs/vcf-representation.md, open point 4).
-- Submitted class as in study_records.sql: the method is null or submitted, and the placement has a position.
-- Binds: study (a study accession).
SELECT
    vs.study_accession,
    vs.variant_study_id,
    vs.variant_accession AS accession,
    p.assembly_accession
FROM dgva.variant_study vs
    JOIN dgva.variant_placement p ON p.variant_study_id = vs.variant_study_id
    JOIN dgva.var_plac_placement_rank r ON r.variant_placement_id = p.variant_placement_id
    LEFT JOIN dgva.variant_placement_method m ON m.variant_placement_id = p.variant_placement_id
WHERE r.placement_rank = 'NotBest'
    AND vs.study_accession = :study
    AND (m.variant_placement_method IS NULL
         OR (m.variant_placement_method NOT LIKE 'Remapped%'
             AND m.variant_placement_method NOT IN ('Remap failed', 'Submitted without placement on assembly',
                                                    'Submitted on non-INSDC assembly')))
    AND NOT EXISTS (SELECT 1 FROM dgva.variant_withdrawn w WHERE w.variant_accession = vs.variant_accession)
