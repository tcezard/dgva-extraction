-- dgva-extraction checks, step "experiments" (ticket 10, point 5): every experiment with what describes it (type,
-- method type, analysis type, detection method, platform), one row per combination, and whether it has calls (one
-- index probe on VARIANT_CALL_EXPERIMENT.EXPERIMENT_ID per experiment). Small tables: 761 experiments. Fast.
-- Binds: none.
SELECT
    e.study_accession,
    e.experiment_id,
    e.submitter_experiment_id,
    e.experiment_type,
    m.method_type,
    a.analysis_type,
    d.detection_method,
    pl.platform_name,
    CASE WHEN EXISTS (SELECT 1 FROM dgva.variant_call_experiment vce WHERE vce.experiment_id = e.experiment_id)
         THEN 1 ELSE 0 END AS has_calls
FROM dgva.dgva_experiment e
    LEFT JOIN dgva.experiment_method em ON em.experiment_id = e.experiment_id
    LEFT JOIN dgva.dgva_method m ON m.method_id = em.method_id
    LEFT JOIN dgva.experiment_analysis ea ON ea.experiment_id = e.experiment_id
    LEFT JOIN dgva.dgva_analysis a ON a.analysis_id = ea.analysis_id
    LEFT JOIN dgva.experiment_detection ed ON ed.experiment_id = e.experiment_id
    LEFT JOIN dgva.dgva_detection d ON d.detection_id = ed.detection_id
    LEFT JOIN dgva.experiment_platform ep ON ep.experiment_id = e.experiment_id
    LEFT JOIN dgva.dgva_platform pl ON pl.platform_id = ep.platform_id
ORDER BY e.study_accession, e.experiment_id
