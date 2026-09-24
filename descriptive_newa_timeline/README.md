# NEWA timeline

Source: Cornell Network for Environment and Weather Applications, [Grape Diseases Model](https://newa.cornell.edu/grape-diseases), Geneva (AgriTech Crittenden), NY (42.8794, -77.0154). Historical tables were retrieved September 21, 2026 and are preserved in `source_logs/`.

`make_timeline.py` selects rows labeled Infection in the Black Rot column, computes grouped-period overlap with each trial window, and creates Figure 6 plus derived CSVs. A period counts when its end is on/after the first spray date and its start precedes the day after final rating. Counts: 21, 10, 10, 14, 13 (2021–2025). No network access is needed.

Grouped periods can include interrupted wetness. They are not confirmed fruit infections. Windows differ among years; actual bloom dates are unavailable. The common-window CSV reports a separate descriptive count of starts, not Figure 6 trial-window totals. Run in a working copy because exports are written beside the script.
