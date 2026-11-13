import { cp, mkdir, rm } from 'node:fs/promises';

const root = new URL('../', import.meta.url);
const webDir = new URL('../www/', import.meta.url);

await rm(webDir, { recursive: true, force: true });
await mkdir(webDir, { recursive: true });

for (const file of ['index.html', 'icon.svg', 'manifest.webmanifest']) {
  await cp(new URL(file, root), new URL(file, webDir));
}

console.log('Prepared local Capacitor web bundle in ./www');
