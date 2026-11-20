-- dgva-extraction checks, step "validations" (ticket 10, point 5): every validation result of the calls validated
-- more than once (about 850,000 of the 2,992,979 rows), to find calls validated by different methods or with mixed
-- results. Full scans of VARIANT_CALL_VALIDATION (3 M rows): run outside busy hours.
-- Binds: none.
SELECT v.variant_call_accession AS accession, v.validation_experiment_id AS experiment_id, v.validation_result
FROM dgva.variant_call_validation v
WHERE v.variant_call_accession IN (
    SELECT variant_call_accession FROM dgva.variant_call_validation
    GROUP BY variant_call_accession HAVING COUNT(*) > 1)
ORDER BY v.variant_call_accession, v.validation_experiment_id
