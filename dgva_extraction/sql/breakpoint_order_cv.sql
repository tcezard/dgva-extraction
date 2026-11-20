-- The values VAR_PLAC_BREAKPOINT_ORDER can take, with their description (ticket 10, open point 4).
-- Binds: none.
SELECT breakpoint_order, description, is_current
FROM dgva.var_plac_breakpoint_order_cv
ORDER BY breakpoint_order
