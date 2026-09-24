#!/usr/bin/env python3
"""Create publication tables from completed sensitivity CSV/JSON outputs.

Only Python's standard library is required. No models or P values are estimated
here: outputs are joined/rescaled and simulation counts are independently checked.
"""
import argparse
import csv
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('base', nargs='?', default='.')
parser.add_argument('--incidence-dir', default='bootstrap_incidence_1999')
parser.add_argument('--severity-dir', default='bootstrap_severity_9999')
args = parser.parse_args()
base = Path(args.base)

def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))

def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

estimate_rows = []
for endpoint, filename, estimate, lower, upper, scale in [
    ('Incidence', 'incidence_EMM_comparison.csv', 'prob', 'asymp.LCL', 'asymp.UCL', 100),
    ('Severity', 'severity_EMM_comparison.csv', 'Estimate_Pct', 'Lower_Pct', 'Upper_Pct', 1),
]:
    source = {r['Program_Class']: r for r in read_csv(base / 'results' / filename)}
    for group in ['Untreated', 'Biological only', 'Mixed', 'Conventional only']:
        r = source[group]
        row = {'Endpoint': endpoint, 'Program_Class': group}
        for dataset in ['Original264', 'Exclude8']:
            row['Estimate_Pct_' + dataset] = scale * float(r[estimate + '_' + dataset])
            row['Lower95_Pct_' + dataset] = scale * float(r[lower + '_' + dataset])
            row['Upper95_Pct_' + dataset] = scale * float(r[upper + '_' + dataset])
        row['Change_Percentage_Points'] = row['Estimate_Pct_Exclude8'] - row['Estimate_Pct_Original264']
        estimate_rows.append(row)
write_csv(base / 'pooled_estimates_table.csv', estimate_rows)

calibration_rows = []
for endpoint, folder, stem in [
    ('Incidence', args.incidence_dir, 'incidence_interaction_score_refit_bootstrap'),
    ('Severity', args.severity_dir, 'severity_interaction_lrt_refit_bootstrap'),
]:
    summary = json.loads((base / folder / (stem + '.json')).read_text())
    original = json.loads((base.parent.parent / 'results' / 'primary_bootstrap' / (stem + '.json')).read_text())
    simulations = read_csv(base / folder / (stem + '_statistics.csv'))
    valid = [float(r['Statistic']) for r in simulations if r['Statistic'] != 'NA']
    assert len(simulations) == summary['requested_simulations']
    assert len(valid) == summary['valid_simulations']
    assert sum(value >= summary['observed_statistic'] for value in valid) == summary['exceedances']
    assert len(simulations) - len(valid) == summary['failed_simulations']
    assert math.isclose(summary['bootstrap_p'], (summary['exceedances'] + 1) / (len(valid) + 1), rel_tol=1e-14)
    assert math.isclose(summary['missing_result_p_lower'], (summary['exceedances'] + 1) / (len(simulations) + 1), rel_tol=1e-14)
    assert math.isclose(summary['missing_result_p_upper'], (summary['exceedances'] + summary['failed_simulations'] + 1) / (len(simulations) + 1), rel_tol=1e-14)
    row = {
        'Endpoint': endpoint,
        'Test': summary['test'],
        'Observed_Statistic_Original264': original['observed_statistic'],
        'Bootstrap_P_Original264': original['bootstrap_p'],
        'Observed_Statistic_Exclude8': summary['observed_statistic'],
        'Requested_Simulations': summary['requested_simulations'],
        'Valid_Simulations': summary['valid_simulations'],
        'Failed_Simulations': summary['failed_simulations'],
        'Exceedances': summary['exceedances'],
        'Bootstrap_P_Exclude8': summary['bootstrap_p'],
        'Monte_Carlo_SE': summary['monte_carlo_se'],
        'Failure_Bound_Lower': summary['missing_result_p_lower'],
        'Failure_Bound_Upper': summary['missing_result_p_upper'],
        'Seed': summary['seed'],
    }
    calibration_rows.append(row)
write_csv(base / 'interaction_calibration_table.csv', calibration_rows)
print('Wrote 8 pooled-estimate rows and 2 interaction-calibration rows; all simulation counts and plus-one P values verified.')
