import SwiftUI
import CoreMotion
import CoreHaptics
import UIKit
import QuartzCore

@MainActor
final class FeedbackModel: ObservableObject {
    @Published var motionOn = false
    @Published var hapticsAvailable = false
    @Published var gravity = CGVector.zero
    @Published var motionHz = 0.0
    @Published var fps = 0.0
    @Published var frameP95 = 0.0
    @Published var hapticDispatchMs = 0.0
    @Published var intensity: Double = 0.65
    @Published var sharpness: Double = 0.45
    @Published var continuousOn = false
    @Published var status = "READY"
    @Published var masses: [CGPoint] = [.zero]

    private let motionManager = CMMotionManager()
    private var engine: CHHapticEngine?
    private var continuousPlayer: CHHapticAdvancedPatternPlayer?
    private var motionTimes: [TimeInterval] = []
    private var frameTimes: [Double] = []
    private var previousFrame: CFTimeInterval = 0
    private var gravityFiltered = CGVector.zero
    private var velocities: [CGVector] = [.zero]
    private var lastSensorWindow: TimeInterval = 0
    private var sensorCount = 0
    private var displayLink: CADisplayLink?
    private var displayTarget: DisplayLinkTarget?

    init() {
        hapticsAvailable = CHHapticEngine.capabilitiesForHardware().supportsHaptics
        prepareHaptics()
    }

    func startMotion() {
        guard motionManager.isDeviceMotionAvailable else { status = "DEVICE MOTION UNAVAILABLE"; return }
        guard !motionOn else { return }
        motionManager.deviceMotionUpdateInterval = 1.0 / 60.0
        motionManager.startDeviceMotionUpdates(using: .xArbitraryCorrectedZVertical, to: .main) { [weak self] data, error in
            guard let self, let data, error == nil else { return }
            Task { @MainActor in self.consumeMotion(data) }
        }
        motionOn = true
        status = "MOTION ACTIVE"
    }

    func stopMotion() {
        motionManager.stopDeviceMotionUpdates()
        motionOn = false
        gravityFiltered = .zero
        gravity = .zero
        status = "MOTION STOPPED"
    }

    private func consumeMotion(_ data: CMDeviceMotion) {
        let g = data.gravity
        let angle = currentInterfaceAngle()
        let c = cos(angle), s = sin(angle)
        // Core Motion gravity is in device coordinates. Project it into the current UI axes.
        let sx = g.x * c + g.y * s
        let sy = -g.x * s + g.y * c
        let dz: (Double) -> Double = { value in
            let a = abs(value), dead = 0.035
            return a < dead ? 0 : (value < 0 ? -1 : 1) * min(1, (a - dead) / (1 - dead))
        }
        let tx = dz(sx), ty = dz(sy)
        gravityFiltered.dx += (tx - gravityFiltered.dx) * 0.16
        gravityFiltered.dy += (ty - gravityFiltered.dy) * 0.16
        gravity = gravityFiltered

        let t = data.timestamp
        sensorCount += 1
        if lastSensorWindow == 0 { lastSensorWindow = t }
        if t - lastSensorWindow >= 1 {
            motionHz = Double(sensorCount) / (t - lastSensorWindow)
            sensorCount = 0
            lastSensorWindow = t
        }
        integrateMasses(dt: motionManager.deviceMotionUpdateInterval)
    }

    private func currentInterfaceAngle() -> Double {
        let scene = UIApplication.shared.connectedScenes
            .compactMap { $0 as? UIWindowScene }.first
        switch scene?.interfaceOrientation {
        case .landscapeLeft: return .pi / 2
        case .landscapeRight: return -.pi / 2
        case .portraitUpsideDown: return .pi
        default: return 0
        }
    }

    private func integrateMasses(dt: TimeInterval) {
        let d = min(0.04, max(0, dt))
        guard d > 0 else { return }
        if velocities.count != masses.count { velocities = Array(repeating: .zero, count: masses.count) }
        for i in masses.indices {
            velocities[i].dx = (velocities[i].dx + gravityFiltered.dx * 190 * d) * exp(-2.4 * d)
            velocities[i].dy = (velocities[i].dy + gravityFiltered.dy * 190 * d) * exp(-2.4 * d)
            masses[i].x += velocities[i].dx * d
            masses[i].y += velocities[i].dy * d
        }
    }

    func addMass() {
        let n = masses.count
        let a = Double(n) * 2.399963
        let r = n == 0 ? 0.0 : min(100.0, 22.0 + Double(n) * 8.0)
        masses.append(CGPoint(x: cos(a) * r, y: sin(a) * r))
        velocities.append(.zero)
        pulse()
    }

    func reset() {
        masses = [.zero]
        velocities = [.zero]
        gravityFiltered = .zero
        gravity = .zero
        stopContinuous()
    }

    private func prepareHaptics() {
        guard hapticsAvailable else { status = "HAPTICS NOT AVAILABLE"; return }
        do {
            let h = try CHHapticEngine()
            h.stoppedHandler = { [weak self] _ in
                Task { @MainActor in self?.status = "HAPTIC ENGINE STOPPED" }
            }
            h.resetHandler = { [weak self] in
                Task { @MainActor in
                    do { try self?.engine?.start() }
                    catch { self?.status = "HAPTIC RESTART FAILED" }
                }
            }
            try h.start()
            engine = h
            status = "CORE HAPTICS READY"
        } catch {
            status = "CORE HAPTICS ERROR: \(error.localizedDescription)"
        }
    }

    func pulse() {
        guard let engine, hapticsAvailable else { status = "HAPTICS UNSUPPORTED"; return }
        let start = CACurrentMediaTime()
        do {
            let events = [
                CHHapticEvent(eventType: .hapticTransient, parameters: [
                    .init(parameterID: .hapticIntensity, value: Float(intensity)),
                    .init(parameterID: .hapticSharpness, value: Float(sharpness))
                ], relativeTime: 0),
                CHHapticEvent(eventType: .hapticTransient, parameters: [
                    .init(parameterID: .hapticIntensity, value: Float(intensity * 0.55)),
                    .init(parameterID: .hapticSharpness, value: Float(sharpness * 0.8))
                ], relativeTime: 0.12)
            ]
            let pattern = try CHHapticPattern(events: events, parameters: [])
            let player = try engine.makePlayer(with: pattern)
            try player.start(atTime: CHHapticTimeImmediate)
            hapticDispatchMs = (CACurrentMediaTime() - start) * 1000
            status = "PATTERN DISPATCHED"
        } catch {
            status = "HAPTIC PLAY FAILED: \(error.localizedDescription)"
        }
    }

    func startContinuous() {
        guard let engine, hapticsAvailable, !continuousOn else { return }
        do {
            let event = CHHapticEvent(eventType: .hapticContinuous, parameters: [
                .init(parameterID: .hapticIntensity, value: Float(intensity)),
                .init(parameterID: .hapticSharpness, value: Float(sharpness))
            ], relativeTime: 0, duration: 60)
            let pattern = try CHHapticPattern(events: [event], parameters: [])
            let player = try engine.makeAdvancedPlayer(with: pattern)
            try player.start(atTime: CHHapticTimeImmediate)
            continuousPlayer = player
            continuousOn = true
            status = "CONTINUOUS HAPTICS ON"
        } catch {
            status = "CONTINUOUS HAPTICS FAILED: \(error.localizedDescription)"
        }
    }

    func updateContinuous() {
        guard continuousOn, let player = continuousPlayer else { return }
        do {
            let parameters = [
                CHHapticDynamicParameter(parameterID: .hapticIntensityControl, value: Float(intensity), relativeTime: 0),
                CHHapticDynamicParameter(parameterID: .hapticSharpnessControl, value: Float(sharpness), relativeTime: 0)
            ]
            try player.sendParameters(parameters, atTime: CHHapticTimeImmediate)
        } catch { status = "HAPTIC PARAMETER UPDATE FAILED" }
    }

    func stopContinuous() {
        guard let player = continuousPlayer else { return }
        do { try player.stop(atTime: CHHapticTimeImmediate) } catch { }
        continuousPlayer = nil
        continuousOn = false
        status = "CONTINUOUS HAPTICS OFF"
    }

    func startFrameMeasurement() {
        guard displayLink == nil else { return }
        let target = DisplayLinkTarget(model: self)
        displayTarget = target
        let link = CADisplayLink(target: target, selector: #selector(DisplayLinkTarget.tick(_:)))
        if #available(iOS 15.0, *) {
            link.preferredFrameRateRange = CAFrameRateRange(minimum: 30, maximum: 120, preferred: 120)
        }
        link.add(to: .main, forMode: .common)
        displayLink = link
    }

    fileprivate func recordFrame(_ timestamp: CFTimeInterval) {
        if previousFrame > 0 {
            frameTimes.append((timestamp - previousFrame) * 1000)
            if frameTimes.count > 1200 { frameTimes.removeFirst(frameTimes.count - 1200) }
            let sorted = frameTimes.sorted()
            fps = sorted.isEmpty ? 0 : 1000 / (sorted.reduce(0, +) / Double(sorted.count))
            frameP95 = sorted.isEmpty ? 0 : sorted[min(sorted.count - 1, Int(Double(sorted.count - 1) * 0.95))]
        }
        previousFrame = timestamp
    }

    func stopAll() {
        stopMotion()
        stopContinuous()
        displayLink?.invalidate()
        displayLink = nil
        displayTarget = nil
        engine?.stop(completionHandler: nil)
    }
}

private final class DisplayLinkTarget: NSObject {
    weak var model: FeedbackModel?
    init(model: FeedbackModel) { self.model = model }
    @objc func tick(_ link: CADisplayLink) { model?.recordFrame(link.timestamp) }
}

struct ContentView: View {
    @StateObject private var model = FeedbackModel()
    var body: some View {
        VStack(spacing: 18) {
            HStack {
                Text("SPACETIME").font(.system(size: 15, weight: .bold, design: .rounded))
                Spacer()
                Text("NATIVE FEEDBACK LAB").font(.system(size: 10, weight: .medium, design: .monospaced)).foregroundStyle(.secondary)
            }
            Canvas { context, size in
                let spacing: CGFloat = 34
                var grid = Path()
                stride(from: CGFloat(0), through: size.width, by: spacing).forEach { x in
                    grid.move(to: CGPoint(x: x, y: 0)); grid.addLine(to: CGPoint(x: x, y: size.height))
                }
                stride(from: CGFloat(0), through: size.height, by: spacing).forEach { y in
                    grid.move(to: CGPoint(x: 0, y: y)); grid.addLine(to: CGPoint(x: size.width, y: y))
                }
                context.stroke(grid, with: .color(.white.opacity(0.38)), lineWidth: 1)
                let center = CGPoint(x: size.width / 2, y: size.height / 2)
                for p in model.masses {
                    let point = CGPoint(x: center.x + p.x, y: center.y + p.y)
                    let rect = CGRect(x: point.x - 5, y: point.y - 5, width: 10, height: 10)
                    context.fill(Path(ellipseIn: rect), with: .color(.cyan.opacity(0.9)))
                    let glow = CGRect(x: point.x - 60, y: point.y - 60, width: 120, height: 120)
                    context.fill(Path(ellipseIn: glow), with: .radialGradient(Gradient(colors: [.cyan.opacity(0.16), .clear]), center: point, startRadius: 2, endRadius: 60))
                }
            }
            .frame(minHeight: 180)
            .background(Color.black)
            .clipShape(RoundedRectangle(cornerRadius: 16))
            .overlay(RoundedRectangle(cornerRadius: 16).stroke(.white.opacity(0.12)))
            HStack(spacing: 9) {
                Button(model.motionOn ? "STOP MOTION" : "ENABLE GRAVITY") {
                    model.motionOn ? model.stopMotion() : model.startMotion()
                }.buttonStyle(.borderedProminent)
                Button("ADD MASS") { model.addMass() }.buttonStyle(.bordered)
                Button("RESET") { model.reset() }.buttonStyle(.bordered)
            }
            HStack(spacing: 9) {
                Button("HAPTIC PATTERN") { model.pulse() }.buttonStyle(.borderedProminent)
                Button(model.continuousOn ? "STOP CONTINUOUS" : "CONTINUOUS") {
                    model.continuousOn ? model.stopContinuous() : model.startContinuous()
                }.buttonStyle(.bordered)
            }
            VStack(alignment: .leading, spacing: 8) {
                HStack { Text("Intensity"); Slider(value: $model.intensity).onChange(of: model.intensity) { _, _ in model.updateContinuous() }; Text(model.intensity.formatted(.number.precision(.fractionLength(2)))) }
                HStack { Text("Sharpness"); Slider(value: $model.sharpness).onChange(of: model.sharpness) { _, _ in model.updateContinuous() }; Text(model.sharpness.formatted(.number.precision(.fractionLength(2)))) }
            }.font(.caption)
            Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 8) {
                GridRow { Text("Motion rate"); Text(String(format: "%.1f Hz", model.motionHz)); Text("Display"); Text(String(format: "%.1f fps", model.fps)) }
                GridRow { Text("Frame p95"); Text(String(format: "%.2f ms", model.frameP95)); Text("Dispatch"); Text(String(format: "%.3f ms", model.hapticDispatchMs)) }
                GridRow { Text("Gravity X"); Text(String(format: "%.3f", model.gravity.dx)); Text("Gravity Y"); Text(String(format: "%.3f", model.gravity.dy)) }
                GridRow { Text("Core Haptics"); Text(model.hapticsAvailable ? "Available" : "Unavailable"); Text("Masses"); Text("\(model.masses.count)") }
            }.font(.system(size: 11, design: .monospaced)).frame(maxWidth: .infinity, alignment: .leading)
            Text(model.status).font(.system(size: 10, weight: .medium, design: .monospaced)).foregroundStyle(.secondary).frame(maxWidth: .infinity, alignment: .leading)
            Spacer(minLength: 0)
        }
        .padding()
        .background(Color(red: 0.025, green: 0.035, blue: 0.055).ignoresSafeArea())
        .preferredColorScheme(.dark)
        .onAppear { model.startFrameMeasurement() }
        .onDisappear { model.stopAll() }
    }
}

@main
struct SpacetimeFeedbackApp: App {
    var body: some Scene { WindowGroup { ContentView() } }
}