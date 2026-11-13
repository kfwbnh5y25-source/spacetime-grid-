# Feedback Lab — native iOS experiment

This Capacitor project is deliberately isolated inside `feedback-lab/`. It does not package, edit, or replace the normal Spacetime app at the repository root.

## Requirements

- macOS with Xcode 26 or newer
- Node.js 22 or newer
- iPhone connected to the Mac, or an iOS Simulator

## First-time setup

From this directory:

```sh
cd feedback-lab
npm install
npm run ios:add
npm run ios:open
```

In Xcode, select the **App** target, set your Apple Development Team under **Signing & Capabilities**, select your iPhone 14 Pro, and press Run.

The native target supports portrait and both landscape orientations. UIKit controls the app-window transition; the Feedback Lab canvas is hosted in WKWebView. No custom rotation animation is injected initially, so we can measure the native system transition before considering a custom native change.

## Updating after Feedback Lab web changes

```sh
cd feedback-lab
npm run ios:sync
npm run ios:open
```

The build script copies only `feedback-lab/index.html`, its service worker and manifest into the native bundle. It rewrites the icon path to a local bundled asset. It never copies the root `index.html`.

## Expected behaviour

- iOS owns orientation rotation and animates the native app window.
- WKWebView receives the native size transition.
- The lab's own JavaScript can record viewport and orientation timing.
- This aims to use the same system-level rotation path as other native apps; exact visual parity with Camera still depends on the content layout and iOS compositor.
