#!/usr/bin/env node
// Usage: node render_publication_figures.mjs FIGURES_DIRECTORY
// Requires Node.js and sharp (npm install sharp). An optional SHARP_MODULE
// environment variable may identify an existing sharp installation.
import fs from 'node:fs/promises';
import path from 'node:path';
import { createRequire } from 'node:module';
const require = createRequire(import.meta.url);
const sharp = require(process.env.SHARP_MODULE || 'sharp');
const dir = path.resolve(process.argv[2] || 'figures');
const names = ['Figure_1_analysis_workflow', 'Figure_2_year_by_program_class',
  'Figure_3_metal_vs_nonmetal', 'Figure_S1_program_level_severity'];
const audit = {};
for (const stem of names) {
  const svg = await fs.readFile(path.join(dir, `${stem}.svg`), 'utf8');
  const viewBox = svg.match(/viewBox="0 0 ([\d.]+) ([\d.]+)"/);
  if (!viewBox) throw new Error(`Missing viewBox for ${stem}`);
  const width = 4200;
  const height = Math.round(width * Number(viewBox[2]) / Number(viewBox[1]));
  // Supply an explicit raster viewport to render native vector geometry directly
  // at 600 dpi and 7in width, avoiding libvips' inch-unit density scaling.
  const viewport = svg.replace(/(<svg\b[^>]*?)width="[^"]+" height="[^"]+"/,
    `$1width="${width}px" height="${height}px"`);
  const pixels = await sharp(Buffer.from(viewport), { density: 72 })
    .flatten({ background: '#FFFFFF' }).removeAlpha().toColourspace('srgb')
    .withMetadata({ density: 600 }).png().toBuffer();
  const meta = await sharp(pixels).metadata();
  if (meta.width !== width || meta.height !== height) throw new Error('Unexpected raster size');
  await fs.writeFile(path.join(dir, `${stem}.png`), pixels);
  await sharp(pixels).withMetadata({ density: 600 }).tiff({ compression: 'lzw' })
    .toFile(path.join(dir, `${stem}.tiff`));
  audit[stem] = {width, height, dpi: 600, width_in: 7, height_in: height/600,
    mode: 'RGB', tiff_compression: 'LZW', raster_upsampling: false};
}
audit.renderer = `sharp ${sharp.versions.sharp}, libvips ${sharp.versions.vips}`;
await fs.writeFile(path.join(dir, 'publication_raster_audit.json'), JSON.stringify(audit, null, 2));
console.log(JSON.stringify(audit, null, 2));
