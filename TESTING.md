Spacetime Grid — device test plan

IPHONE 14 PRO (SAFARI)
1. Open the site in Safari, then Share → Add to Home Screen. Launch from the new icon. Expected: standalone window without Safari address/title bars; black background fills the screen and content respects Dynamic Island and home indicator.
2. Swipe with one finger. Expected: existing grid deformation still works.
3. Hold one finger still for about 0.6 s. Expected: a persistent mass is placed and a brief MASS PLACED label appears. It remains after reload.
4. Use two fingers: pinch in/out and drag. Expected: zoom (0.7x–2.5x) and pan; release saves the view.
5. Reload the app. Expected: camera view and placed masses restore.
6. Add more than eight masses. Expected: oldest mass is replaced.
7. Swipe the app away, reopen it. Expected: animation resumes; no runaway background rendering.
8. Enable iOS Settings → Accessibility → Motion → Reduce Motion. Expected: animated wave phase stops; touch deformation remains.
9. Load once online, enable Airplane Mode, relaunch. Expected: app shell loads offline; remote starfield may be absent on a cold cache.
10. Rotate portrait ↔ landscape. Expected: grid resizes without clipping or a stuck transform.

MAC (SAFARI)
1. Open the site; move the pointer. Expected: planet cursor follows smoothly and grid bends.
2. Click anywhere. Expected: no click-created deformation/action.
3. Press Space. Expected: wave animation pauses; press Space again to resume.
4. Press + / - to zoom, arrow keys to pan, R to reset view and clear placed masses.
5. Press S. Expected: share sheet when supported, otherwise copy/share fallback.
6. Use browser zoom and resize the window. Expected: canvas remains full-window.
7. Add to Dock from Safari if desired; launch standalone and verify app-like window.

PERFORMANCE / REGRESSION
1. Run the benchmark page at query ?bench=1 in the foreground for 60 s at DPR 2 and cap 60; repeat while interacting.
2. Confirm WebGL 2 renderer, stable frame cadence, no long tasks, and no unexpected CPU increase.
3. Repeat on iPhone with DPR 3 and 60 FPS; compare temperature and frame intervals to the previous baseline.
4. Test offline launch after one successful online load.
5. Verify existing visual quality: no bright segment junctions, no dotted/dashed grid lines, no iOS text loupe, and no Mac click interaction.

PLATFORM LIMITATIONS
- iOS Safari does not expose a general-purpose web API for the native Taptic Engine; reliable custom haptics are not promised.
- Web pages cannot force the display refresh rate or continue animation freely in the background. The page pauses rendering when hidden.
- iOS Home Screen icon support for SVG varies by iOS version; the manifest icon is SVG. If iOS ignores it, Safari may generate an icon from the page.
