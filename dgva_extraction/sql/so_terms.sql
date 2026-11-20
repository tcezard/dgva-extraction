-- Ticket 10: the variant types in use and their SO terms, to check that etc/symbolic_alleles.yaml has an entry for
-- each (dgva-extraction check-types). Counts are over all variants, withdrawn or not.
-- Binds: none.
SELECT 'region' AS region_or_call, t.exchange_variant_region_type AS exchange_type, t.sequence_ontology_term AS so_term,
    so.sequence_ontology_id AS so_id, COUNT(r.variant_region_accession) AS variants
FROM dgva.variant_region_type_cv t
    LEFT JOIN dgva.sequence_ontology_cv so ON so.sequence_ontology_term = t.sequence_ontology_term
    LEFT JOIN dgva.variant_region r ON r.exchange_variant_region_type = t.exchange_variant_region_type
GROUP BY t.exchange_variant_region_type, t.sequence_ontology_term, so.sequence_ontology_id
UNION ALL
SELECT 'call', t.exchange_variant_call_type, t.sequence_ontology_term, so.sequence_ontology_id,
    COUNT(c.variant_call_accession)
FROM dgva.variant_call_type_cv t
    LEFT JOIN dgva.sequence_ontology_cv so ON so.sequence_ontology_term = t.sequence_ontology_term
    LEFT JOIN dgva.variant_call c ON c.exchange_variant_call_type = t.exchange_variant_call_type
GROUP BY t.exchange_variant_call_type, t.sequence_ontology_term, so.sequence_ontology_id
ORDER BY 1, 2
