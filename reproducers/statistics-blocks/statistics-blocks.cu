// Cost per call of the handling of the statistics of a synchronous batch
// without caller statistics, for several designs. The kernel mimics the
// kernel of MFront/GenericParallel/GPU.hxx: each thread computes the status
// of its point, the statistics are only accessed if a point of the warp
// reports an event.
//
// V0 bare      : kernel + stream synchronization (lower bound)
// V1 copies    : H2D copy from a host local, kernel, sync, D2H copy, sync
//                (TFEL master)
// V2 readback  : device block kept neutral, kernel, async D2H copy to pinned
//                memory, sync
// V3 flag      : device block kept neutral, flag in pinned memory written by
//                the threads having an event, sync, flag read by the host
// V4 managed   : block in managed memory initialized and read by the host
//                (former prototype)
//
// Each variant is measured alone, then after a kernel of ~200 us enqueued on
// the same stream before each call (work of the calling code, for instance
// the computation of the gradients), which shows the hidden synchronizations.
#include <cstdio>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cuda_runtime.h>

#define CHECK(x)                                                        \
  do {                                                                  \
    const auto e = (x);                                                 \
    if (e != cudaSuccess) {                                             \
      std::printf("%s:%d %s\n", __FILE__, __LINE__, cudaGetErrorString(e)); \
      std::exit(1);                                                     \
    }                                                                   \
  } while (0)

struct Statistics {
  unsigned long long number_of_failures;
  unsigned long long first_failure;
  unsigned long long number_of_unreliable_results;
  unsigned long long first_unreliable_result;
  unsigned long long rdt;
};

// event at point i if fail_every != 0 and i % fail_every == 0
__global__ void apply(const double* const x, double* const y,
                      const unsigned long long n, Statistics* const s,
                      int* const flag, const unsigned long long fail_every) {
  const auto i = static_cast<unsigned long long>(blockIdx.x) * blockDim.x +
                 threadIdx.x;
  bool failed = false;
  if (i < n) {
    y[i] = 2 * x[i] + 1;
    failed = (fail_every != 0) && (i % fail_every == 0);
  }
  if (!__any_sync(0xffffffff, failed)) {
    return;
  }
  const auto nf = __popc(__ballot_sync(0xffffffff, failed));
  if (threadIdx.x % 32 == 0) {
    atomicAdd(&(s->number_of_failures), static_cast<unsigned long long>(nf));
    atomicMin(&(s->first_failure), i);
    if (flag != nullptr) {
      *reinterpret_cast<volatile int*>(flag) = 1;
    }
  }
}

__global__ void busy(double* const y, const int n, const long long cycles) {
  const auto start = clock64();
  while (clock64() - start < cycles) {
  }
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    y[0] += 1;
  }
}

static Statistics neutral() {
  return {0, ~0ull, 0, ~0ull, ~0ull};
}

int main(int argc, char** argv) {
  const int repetitions = 2000;
  cudaStream_t stream;
  CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
  const unsigned long long nmax = 1u << 22;
  double *x, *y;
  CHECK(cudaMalloc(&x, nmax * sizeof(double)));
  CHECK(cudaMalloc(&y, nmax * sizeof(double)));
  CHECK(cudaMemset(x, 0, nmax * sizeof(double)));
  Statistics* d_block;
  CHECK(cudaMalloc(&d_block, sizeof(Statistics)));
  Statistics* m_block;
  CHECK(cudaMallocManaged(&m_block, sizeof(Statistics)));
  Statistics* p_readback;
  CHECK(cudaMallocHost(&p_readback, sizeof(Statistics)));
  int* p_flag;
  CHECK(cudaMallocHost(&p_flag, sizeof(int)));
  *p_flag = 0;
  {
    const auto s = neutral();
    CHECK(cudaMemcpy(d_block, &s, sizeof(Statistics), cudaMemcpyHostToDevice));
  }
  // ~200 us of work of the calling code
  int clock_rate = 0;
  CHECK(cudaDeviceGetAttribute(&clock_rate, cudaDevAttrClockRate, 0));
  const long long busy_cycles = static_cast<long long>(clock_rate) * 200 / 1000;

  auto launch = [&](unsigned long long n, Statistics* s, int* f,
                    unsigned long long fe) {
    const unsigned threads = 128;
    const auto blocks = static_cast<unsigned>((n + threads - 1) / threads);
    apply<<<blocks, threads, 0, stream>>>(x, y, n, s, f, fe);
  };

  auto variant = [&](int v, unsigned long long n, unsigned long long fe,
                     bool pipeline) {
    unsigned long long failures = 0;
    if (pipeline) {
      busy<<<1, 32, 0, stream>>>(y, 1, busy_cycles);
    }
    if (v == 0) {
      launch(n, d_block, nullptr, fe);
      CHECK(cudaStreamSynchronize(stream));
    } else if (v == 1) {
      auto h = neutral();
      CHECK(cudaMemcpyAsync(d_block, &h, sizeof(Statistics),
                            cudaMemcpyHostToDevice, stream));
      launch(n, d_block, nullptr, fe);
      CHECK(cudaStreamSynchronize(stream));
      CHECK(cudaMemcpyAsync(&h, d_block, sizeof(Statistics),
                            cudaMemcpyDeviceToHost, stream));
      CHECK(cudaStreamSynchronize(stream));
      failures = h.number_of_failures;
    } else if (v == 2) {
      launch(n, d_block, nullptr, fe);
      CHECK(cudaMemcpyAsync(p_readback, d_block, sizeof(Statistics),
                            cudaMemcpyDeviceToHost, stream));
      CHECK(cudaStreamSynchronize(stream));
      failures = p_readback->number_of_failures;
      if (failures != 0) {
        const auto s = neutral();
        CHECK(cudaMemcpy(d_block, &s, sizeof(Statistics),
                         cudaMemcpyHostToDevice));
      }
    } else if (v == 3) {
      launch(n, d_block, p_flag, fe);
      CHECK(cudaStreamSynchronize(stream));
      if (*reinterpret_cast<volatile int*>(p_flag) != 0) {
        Statistics h;
        CHECK(cudaMemcpy(&h, d_block, sizeof(Statistics),
                         cudaMemcpyDeviceToHost));
        failures = h.number_of_failures;
        const auto s = neutral();
        CHECK(cudaMemcpy(d_block, &s, sizeof(Statistics),
                         cudaMemcpyHostToDevice));
        *p_flag = 0;
      }
    } else if (v == 4) {
      *m_block = neutral();
      launch(n, m_block, nullptr, fe);
      CHECK(cudaStreamSynchronize(stream));
      failures = m_block->number_of_failures;
    }
    return failures;
  };

  const char* names[] = {"V0 bare", "V1 copies", "V2 readback", "V3 flag",
                         "V4 managed"};
  for (const auto pipeline : {false, true}) {
    std::printf("\n=== %s\n", pipeline ? "after 200 us of work on the stream"
                                         : "alone");
    std::printf("%-12s", "n");
    for (const auto* nm : names) {
      std::printf("%14s", nm);
    }
    std::printf("   (us per call, median of %d)\n", repetitions);
    for (const auto n : {1ull, 1000ull, 100000ull, 1ull << 22}) {
      std::printf("%-12llu", n);
      for (int v = 0; v != 5; ++v) {
        for (int k = 0; k != 50; ++k) {
          variant(v, n, 0, pipeline);
        }
        std::vector<double> t(repetitions);
        for (int k = 0; k != repetitions; ++k) {
          const auto t0 = std::chrono::steady_clock::now();
          variant(v, n, 0, pipeline);
          const auto t1 = std::chrono::steady_clock::now();
          t[k] = std::chrono::duration<double, std::micro>(t1 - t0).count();
        }
        std::nth_element(t.begin(), t.begin() + repetitions / 2, t.end());
        std::printf("%14.1f", t[repetitions / 2]);
      }
      std::printf("\n");
    }
  }
  // correctness of the designs keeping the block neutral: failures, then
  // success, then failures
  std::printf("\n=== check (expected failures: 8, 0, 8)\n");
  for (int v = 1; v != 5; ++v) {
    const auto a = variant(v, 1000, 125, false);
    const auto b = variant(v, 1000, 0, false);
    const auto c = variant(v, 1000, 125, false);
    std::printf("%-12s %llu %llu %llu\n", names[v], a, b, c);
  }
  return 0;
}
