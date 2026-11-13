# Native iOS build (Capacitor)

The GitHub Pages/PWA build remains unchanged. Capacitor is a separate native target that packages the same `index.html` locally in a WKWebView.

## Requirements

- macOS with Xcode 26 or newer
- Node.js 22 or newer
- iPhone connected to the Mac (or an iOS Simulator)

## First-time setup

From the repository root:

```sh
npm install
npm run ios:add
npm run ios:open
```

In Xcode, select the **App** target, set your Apple Development Team under **Signing & Capabilities**, select your iPhone, then press Run.

The iOS target allows portrait and both landscape orientations. UIKit owns the app-window rotation and WKWebView receives the resulting native layout changes. The app's grid renderer itself is not rotated to fake the device transition.

## Updating the native app after web changes

```sh
npm run ios:sync
npm run ios:open
```

## Notes

- This does not change or replace the GitHub Pages web app.
- The native app uses bundled local assets, not a remote URL, so it does not depend on network availability to launch.
- No custom rotation animation is added initially. First test UIKit's native orientation transition; only add Swift-level intervention if the real device shows a defect.
- iPhone upside-down portrait is intentionally not enabled; iOS apps generally do not support it on iPhone.
- Capacitor uses WKWebView and requires Xcode for the iOS build/signing process.
