#if os(macOS)
import Darwin
#else
import WinSDK
#endif

// Suitability of QueryPerformanceCounter for FP64 numbers:

/*
import WinSDK

var frequency = LARGE_INTEGER()
let boolean = QueryPerformanceCounter(&frequency)

print(frequency.LowPart as UInt32)
print(frequency.HighPart as Int32)
print(frequency.QuadPart as Int64)
print(1 << 53)

// 2085122834
// 17871
// 76757445669650
// 9007199254740992
*/

struct Profiler {
  static func time() -> Double {
    func queryTickCount() -> UInt64 {
      #if os(macOS)
      return mach_continuous_time()
      #else
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
