# Spacetime Feedback Lab

This is a separate experiment. It does not modify the root Spacetime simulation.

## Web app
Open \`https://kfwbnh5y25-source.github.io/spacetime-grid-/feedback-lab/\` in Safari, then use Share → Add to Home Screen for the standalone PWA. The app has no CDN, analytics, libraries, or remote runtime dependencies.

The PWA measures:
- requestAnimationFrame cadence: sample count, mean, p50, p95, max and estimated FPS;
- input-to-next-frame delay (software timing, not display scan-out);
- device-motion sample rate and gravity projected into the active screen axes;
- orientation changes and the screen angle;
- IndexedDB write round-trip;
- API exposure and dispatch timing for vibration, share and wake lock;
- online/offline and service-worker status.

Tap **10s Measure** with the app foregrounded. Export JSON to save the report. These are instrumented live measurements, not pre-filled benchmark claims. Physical haptic amplitude and actual display scan-out latency cannot be measured by JavaScript alone.

### iPhone limitations
- iOS motion access is requested from the **Enable Gravity** button because Safari requires a user gesture.
- The PWA tests \`navigator.vibrate()\` only when present. Current iOS Safari does not expose it; no false “haptic succeeded” result is shown.
- Screen Wake Lock is requested only after the user taps **Keep Awake**, and is reacquired only when the user had it enabled.
- High refresh is measured from the actual animation-frame cadence; the app does not assume 120 Hz merely because the iPhone display supports ProMotion.

## Native iOS companion
\`native/SpacetimeFeedbackApp.swift\` is a dependency-free SwiftUI app source using only Apple frameworks (Core Motion, Core Haptics, SwiftUI, UIKit and QuartzCore).

To run it, create a new iOS App project in Xcode (SwiftUI lifecycle), add the Swift file to the target, and add \`NSMotionUsageDescription\` to Info.plist using the string in \`Info.plist\`. Select an iPhone 14 Pro and run on device. The iOS Simulator cannot validate physical haptics.

The native test reports Core Motion sample rate, display-link cadence, and software haptic dispatch time. Core Haptics supports transient patterns and continuous haptics with adjustable intensity/sharpness on compatible hardware. Software dispatch time is not the same as measured physical actuator onset.

## Test protocol
1. Run the PWA in portrait for 10 seconds; export JSON.
2. Rotate through both landscape orientations while masses are moving. The masses stay in center-relative scene coordinates; only the toolbar changes edge. The gravity vector is reprojected after orientation changes.
3. Repeat the 10-second capture in landscape; compare FPS, p95 frame interval, motion Hz and orientation count.
4. Repeat steps 1–3 in the native companion. Compare measured values, not assumed capabilities.
5. Test offline by loading once online, enabling Airplane Mode, then reopening the app. A previously cached shell should still load.

## Interaction
- Drag one pointer to pan; two pointers pinch, twist and translate.
- On Mac, arrow keys pan smoothly and + / - zoom smoothly. Trackpad pinch/twist uses Safari GestureEvents.
- Add mass creates a persistent test particle. Reset clears camera and particles.
- Enable Gravity requests motion permission; masses experience a damped gravity-like force in the current screen axes.
- Haptics are deliberately separated: browser vibration where exposed, Core Haptics in the native companion.
