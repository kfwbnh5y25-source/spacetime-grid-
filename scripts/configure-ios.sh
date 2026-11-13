#!/bin/bash
set -euo pipefail

PLIST="ios/App/App/Info.plist"
if [[ ! -f "$PLIST" ]]; then
  echo "Missing $PLIST. Run: npm run ios:add"
  exit 1
fi

PB=/usr/libexec/PlistBuddy
"$PB" -c "Delete :UISupportedInterfaceOrientations" "$PLIST" 2>/dev/null || true
"$PB" -c "Add :UISupportedInterfaceOrientations array" "$PLIST"
"$PB" -c "Add :UISupportedInterfaceOrientations:0 string UIInterfaceOrientationPortrait" "$PLIST"
"$PB" -c "Add :UISupportedInterfaceOrientations:1 string UIInterfaceOrientationLandscapeLeft" "$PLIST"
"$PB" -c "Add :UISupportedInterfaceOrientations:2 string UIInterfaceOrientationLandscapeRight" "$PLIST"

echo "Configured iPhone portrait + both landscape orientations."
