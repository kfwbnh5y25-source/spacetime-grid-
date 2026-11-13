import { cp, mkdir, readFile, rm, writeFile } from 'node:fs/promises';

const root = new URL('../', import.meta.url);
const webDir = new URL('../www/', import.meta.url);
await rm(webDir, { recursive: true, force: true });
await mkdir(webDir, { recursive: true });

await cp(new URL('index.html', root), new URL('index.html', webDir));
const serviceWorker = await readFile(new URL('sw.js', root), 'utf8');
await writeFile(new URL('sw.js', webDir), serviceWorker.replaceAll('../icon.svg', './icon.svg'));

const manifest = JSON.parse(await readFile(new URL('manifest.webmanifest', root), 'utf8'));
manifest.icons = [{ src: './icon.svg', sizes: 'any', type: 'image/svg+xml', purpose: 'any maskable' }];
await writeFile(new URL('manifest.webmanifest', webDir), JSON.stringify(manifest, null, 2) + '\n');
await cp(new URL('../../icon.svg', root), new URL('icon.svg', webDir));
console.log('Prepared Feedback Lab only in feedback-lab/www');
