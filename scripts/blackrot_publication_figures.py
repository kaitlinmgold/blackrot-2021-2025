#!/usr/bin/env python3
"""Render publication SVG plots from saved estimates, without fitting any models.

Usage: python3 blackrot_publication_figures.py RESULTS_DIRECTORY OUTPUT_DIRECTORY
Python standard library only. Each SVG has a 7-inch physical width; one viewBox
unit equals one printed point. Rasterize the resulting SVG files at 600 dpi with
render_publication_figures.mjs (Node.js and sharp), or another SVG renderer.
"""
import csv
import hashlib
import html
import json
import re
import sys
from pathlib import Path

CLASS = ['Untreated', 'Biological only', 'Mixed', 'Conventional only']
COLORS = dict(zip(CLASS, ['#4D4D4D', '#009E73', '#0072B2', '#D55E00']))
SHAPES = dict(zip(CLASS, ['circle', 'square', 'triangle', 'diamond']))
METAL = ['Untreated', 'Non-metal biological', 'Metal-containing biological']
MCOLORS = dict(zip(METAL, ['#4D4D4D', '#009E73', '#CC79A7']))
MSHAPES = dict(zip(METAL, ['circle', 'square', 'open-square']))


class SVG:
    def __init__(self, height, title):
        self.height = height
        self.parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="7in" height="{height / 72:g}in" viewBox="0 0 504 {height}">',
                      f'<title>{html.escape(title)}</title>',
                      '<desc>Scientific plot drawn directly from stored model estimates and confidence limits; no statistical models are refitted.</desc>',
                      '<rect width="504" height="100%" fill="white"/>',
                      '<g font-family="Arial, Helvetica, sans-serif" fill="#111111">']
        self.font_sizes = []
        self.points = []

    def text(self, x, y, s, size=8.2, weight='normal', anchor='start', rotate=None):
        self.font_sizes.append(size)
        transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ''
        self.parts.append(f'<text x="{x:.4f}" y="{y:.4f}" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}"{transform}>{html.escape(str(s))}</text>')

    def line(self, x1, y1, x2, y2, color='#606060', width=.6):
        self.parts.append(f'<line x1="{x1:.5f}" y1="{y1:.5f}" x2="{x2:.5f}" y2="{y2:.5f}" stroke="{color}" stroke-width="{width}"/>')

    def marker(self, x, y, shape, color, r=2.35):
        if shape == 'circle':
            element = f'<circle cx="{x:.5f}" cy="{y:.5f}" r="{r}" fill="{color}"/>'
        elif shape in ('square', 'open-square'):
            fill = 'white' if shape == 'open-square' else color
            element = f'<rect x="{x-r:.5f}" y="{y-r:.5f}" width="{2*r}" height="{2*r}" fill="{fill}" stroke="{color}" stroke-width=".9"/>'
        elif shape == 'triangle':
            element = f'<path d="M{x:.5f},{y-r*1.3:.5f} L{x+r*1.25:.5f},{y+r:.5f} L{x-r*1.25:.5f},{y+r:.5f} Z" fill="{color}"/>'
        else:
            element = f'<path d="M{x:.5f},{y-r*1.3:.5f} L{x+r*1.25:.5f},{y:.5f} L{x:.5f},{y+r*1.3:.5f} L{x-r*1.25:.5f},{y:.5f} Z" fill="{color}"/>'
        self.parts.append(element)

    def save(self, target):
        assert min(self.font_sizes) >= 8
        target.write_text('\n'.join(self.parts + ['</g></svg>']), encoding='utf-8')
        return {'physical_inches': [7, self.height / 72], 'minimum_font_pt': min(self.font_sizes),
                'panel_letter_pt': 14 if self.height == 360 else None,
                'points': self.points}


def read_csv(p, used):
    used[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    with p.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def read_p(p, used):
    used[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return float(json.loads(p.read_text())['bootstrap_p'])


def p_phrase(prefix, p):
    return f'{prefix} < 0.001' if p < .001 else f'{prefix} = {p:.3f}'


def percent_rows(rows, is_incidence=False):
    result = []
    for row in rows:
        r = dict(row)
        r['Year'] = int(r['Year'])
        for out, source in [('Estimate_Pct', 'prob'), ('Lower_Pct', 'asymp.LCL'), ('Upper_Pct', 'asymp.UCL')]:
            r[out] = 100 * float(row[source]) if is_incidence else float(row[out])
        assert 0 <= r['Lower_Pct'] <= r['Estimate_Pct'] <= r['Upper_Pct'] <= 100
        result.append(r)
    return result


def top_legend(svg, levels, colors, shapes, positions):
    for x, group in zip(positions, levels):
        svg.line(x-5, 20, x+5, 20, colors[group], .9)
        svg.marker(x, 20, shapes[group], colors[group])
        svg.text(x+8, 22.8, group, size=8.2)


def main_panel(svg, rows, left, title, letter, lines, levels, colors, shapes, group_key, years):
    top, bottom, width = 102, 318, 184
    svg.text(left-42, 54, letter, size=14, weight='bold')
    svg.text(left, 54, title, size=10.5, weight='bold')
    for k, s in enumerate(lines):
        svg.text(left, 72 + k * 13, s, size=8)
    sy = lambda v: bottom - v / 103 * (bottom-top)
    sx = lambda v: left + (v + .6) / (len(years)-1 + 1.2) * width
    for y in range(0, 101, 20):
        svg.line(left, sy(y), left+width, sy(y), '#E6E6E6', .7)
        svg.line(left-2, sy(y), left, sy(y), '#555555')
        svg.text(left-4, sy(y)+2.8, y, anchor='end')
    for i, year in enumerate(years):
        svg.line(sx(i), bottom, sx(i), bottom+3, '#555555')
        svg.text(sx(i), bottom+13, year, anchor='middle')
    for a, b, c, d in [(left,top,left+width,top),(left,bottom,left+width,bottom),(left,top,left,bottom),(left+width,top,left+width,bottom)]:
        svg.line(a,b,c,d,'#707070',.6)
    svg.text(left+width/2, 349, 'Trial year', size=9, anchor='middle')
    ylabel = 'Clusters with symptoms (%)' if title == 'Incidence' else 'Mean diseased berries (%)'
    svg.text(left-34, (top+bottom)/2, ylabel, size=9, anchor='middle', rotate=-90)
    offsets = [.55 * ((i+.5)/len(levels)-.5) for i in range(len(levels))]
    seen = set()
    for row in rows:
        group, year = row[group_key], row['Year']
        assert (year, group) not in seen
        seen.add((year, group))
        g = levels.index(group)
        x = sx(years.index(year) + offsets[g])
        est, low, high = (row[k] for k in ['Estimate_Pct','Lower_Pct','Upper_Pct'])
        svg.line(x, sy(low), x, sy(high), colors[group], .9)
        cap = width / (len(years)-1+1.2) * .08 / len(levels) / 2
        svg.line(x-cap, sy(low), x+cap, sy(low), colors[group], .9)
        svg.line(x-cap, sy(high), x+cap, sy(high), colors[group], .9)
        svg.marker(x, sy(est), shapes[group], colors[group])
        svg.points.append({'endpoint': title, 'year': year, 'group': group,
                           'estimate_pct': est, 'lower_pct': low, 'upper_pct': high})
    assert len(seen) == len(years)*len(levels)


def program_panel(svg, frame, left, top):
    width, height = 128, 197
    year = int(frame[0]['Year'])
    svg.text(left, top-13, f'{year} (n = {len(frame)} programs)', size=9.2, weight='bold')
    sx = lambda v: left + v / 103 * width
    sy = lambda i: top + (i+.5)/len(frame)*height
    for i, row in enumerate(frame):
        y = sy(i)
        svg.line(left, y, left+width, y, '#E8E8E8', .6)
        code = re.sub(r'^Y[0-9]{4}_', '', row['Public_Program_ID'])
        svg.text(left-4, y+2.7, code, size=8, anchor='end')
        svg.line(left-2, y, left, y, '#555555')
        group = row['Program_Class']
        est, lo, hi = [float(row[k]) for k in ['Severity_Mean_Pct','Severity_Lower95_Pct','Severity_Upper95_Pct']]
        assert 0 <= lo <= est <= hi <= 100
        svg.line(sx(lo), y, sx(hi), y, '#8C8C8C', .75)
        svg.marker(sx(est), y, SHAPES[group], COLORS[group], r=2.05)
        svg.points.append({'year': year, 'program': row['Public_Program_ID'], 'group': group,
                           'estimate_pct': est, 'lower_pct': lo, 'upper_pct': hi})
    for y in [top, top+height]: svg.line(left,y,left+width,y,'#707070')
    for x in [left,left+width]: svg.line(x,top,x,top+height,'#707070')
    for v in range(0,101,20):
        svg.line(sx(v),top+height,sx(v),top+height+3,'#555555')
        svg.text(sx(v),top+height+13,v,size=8,anchor='middle')
    svg.text(left+width/2,top+height+28,'Mean cluster severity (%)',size=8.2,anchor='middle')


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    results, out = map(Path, sys.argv[1:])
    out.mkdir(parents=True, exist_ok=True)
    used, audit = {}, {}
    inc = percent_rows(read_csv(results/'primary/incidence_descriptive_year_estimates.csv',used), True)
    sev = percent_rows(read_csv(results/'primary/severity_year_estimates_kenward_roger.csv',used))
    pi = read_p(results/'primary_bootstrap/incidence_interaction_score_refit_bootstrap.json',used)
    ps = read_p(results/'primary_bootstrap/severity_interaction_lrt_refit_bootstrap.json',used)
    pc = read_p(results/'primary_bootstrap/severity_class_lrt_refit_bootstrap.json',used)
    fig = SVG(360, 'Year-specific incidence and severity estimates by operational program class')
    top_legend(fig,CLASS,COLORS,SHAPES,[43,132,248,325])
    main_panel(fig,inc,50,'Incidence','A',['Class: P < 0.001;',p_phrase('year × class: bootstrap P',pi)],CLASS,COLORS,SHAPES,'Program_Class',list(range(2021,2026)))
    main_panel(fig,sev,306,'Severity','B',[p_phrase('Class: bootstrap P',pc)+';',p_phrase('year × class: bootstrap P',ps)],CLASS,COLORS,SHAPES,'Program_Class',list(range(2021,2026)))
    audit['Figure_2_year_by_program_class'] = fig.save(out/'Figure_2_year_by_program_class.svg')
    audit['Figure_2_year_by_program_class']['p_values'] = {'incidence_interaction':pi,'severity_interaction':ps,'severity_class':pc}

    mi = percent_rows(read_csv(results/'secondary/metal_any_incidence_year_estimates.csv',used),True)
    ms = percent_rows(read_csv(results/'secondary/metal_any_severity_year_estimates.csv',used))
    def metal_p(endpoint):
        rows=read_csv(results/f'secondary/metal_any_{endpoint}_pooled_contrasts.csv',used)
        vals=[float(r['p.value']) for r in rows if re.search(r'Non-metal.*Metal-containing',r['contrast'])]
        assert len(vals)==1
        return vals[0]
    mpi,mps=metal_p('incidence'),metal_p('severity')
    fig=SVG(360,'Year-specific incidence and severity estimates in the biological metal comparison')
    top_legend(fig,METAL,MCOLORS,MSHAPES,[40,132,292])
    main_panel(fig,mi,50,'Incidence','A',['Non-metal vs. metal, pooled:',p_phrase('P',mpi)],METAL,MCOLORS,MSHAPES,'Analysis_Group',[2023,2024,2025])
    main_panel(fig,ms,306,'Severity','B',['Non-metal vs. metal, pooled:',p_phrase('P',mps)],METAL,MCOLORS,MSHAPES,'Analysis_Group',[2023,2024,2025])
    audit['Figure_3_metal_vs_nonmetal']=fig.save(out/'Figure_3_metal_vs_nonmetal.svg')
    audit['Figure_3_metal_vs_nonmetal']['p_values']={'incidence_pooled':mpi,'severity_pooled':mps}

    rows=read_csv(results/'descriptive/program_summary.csv',used)
    fig=SVG(552,'Program-level mean cluster severity and confidence limits in five annual trials')
    for year,left,top in [(2021,32,38),(2022,198,38),(2023,364,38),(2024,32,311),(2025,198,311)]:
        frame=[r for r in rows if int(r['Year'])==year]
        program_panel(fig,frame,left,top)
    for j,group in enumerate(CLASS):
        fig.marker(367,390+j*25,SHAPES[group],COLORS[group],2.4)
        fig.text(377,393+j*25,group,size=8.2)
    audit['Figure_S1_program_level_severity']=fig.save(out/'Figure_S1_program_level_severity.svg')
    assert len(audit['Figure_S1_program_level_severity']['points'])==66
    audit['source_sha256']={str(Path(p).relative_to(results)): digest for p,digest in used.items()}
    audit['python_version']=sys.version
    audit['method']='Direct SVG geometry from saved CSVs; percent conversion only for incidence; no model fitting, data exclusion, confidence-limit recomputation, or raster resampling.'
    (out/'publication_plot_audit.json').write_text(json.dumps(audit,indent=2))
    print('Created Figures 2, 3, and S1 SVGs; all font sizes >=8pt at 7in. Numeric inputs recorded in publication_plot_audit.json.')


if __name__=='__main__': main()
