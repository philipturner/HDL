#if os(macOS)
import Darwin
#else
import WinSDK
#endif

struct Profiler {
  static func time() -> Double {
    func queryTickCount() -> UInt64 {
      #if os(macOS)
      return mach_continuous_time()
      #else
      // Hope this is robust on Windows...
      var largeInteger = LARGE_INTEGER()
      QueryPerformanceCounter(&largeInteger)
      return UInt64(largeInteger.QuadPart)
      #endif
    }
    
    func ticksPerSecond() -> Int {
      #if os(macOS)
      return 24_000_000
      #else
      return 10_000_000
      #endif
    }
    
    let elapsedTicks = queryTickCount()
    return Double(elapsedTicks) / Double(ticksPerSecond())
  }
}
